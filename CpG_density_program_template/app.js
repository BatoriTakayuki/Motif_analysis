/* UI controller. Scientific calculations live in core.js / server.py. */
"use strict";
(() => {
  const $ = id => document.getElementById(id);
  const C = window.MotifCore;
  const fmt = n => Number(n).toLocaleString("en-US");
  const number = n => n == null ? "—" : String(Number(Number(n).toPrecision(5)));
  const esc = s => String(s ?? "").replace(/[&<>"']/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
  const state = {source:"sequence", refmode:"paste", mode:"sliding", server:false, busy:false, reference:null, genome:null, location:null, locationCandidates:[], bw:null, seq:null, scan:null, scanMotif:null, guideResult:null, result:null, stale:false, page:0, seqPage:0, chartChrom:null, svg:null};
  const statuses = {ok:"有効", ambiguous_bases:"曖昧塩基あり", chromosome_not_found:"染色体なし", outside_reference:"参照範囲外", no_data:"データなし"};
  const bwDisplayName=()=>$("bw-label").value.trim().replace(/[\r\n\t]/g," ").slice(0,40)||"bigWig";
  const axisLabel=text=>Array.from(text).length>25?Array.from(text).slice(0,24).join("")+"…":text;
  function configText(config){return config.mode==="raw"?"元のbigWig値（追加平滑化なし）":config.mode==="bed"?"BED領域":config.mode==="bin"?`Bin · ${fmt(config.bin_size)} bp（末尾は実長）`:`Sliding window · ${fmt(2*config.flank+1)} bp / ${fmt(config.step)} bp step`;}
  function tableData(result=state.result){
    if(result?.independent)return state.tableTrack==="bigwig"?{rows:result.bwRows,source:"bigwig",mode:result.bwConfig.mode}:{rows:result.rows,source:"sequence",mode:result.mode};
    return result;
  }
  function refreshFigureLabels(result){
    const name=bwDisplayName(),overlay=result.source==="overlay";
    $("chart-title").textContent=overlay?`${result.motif}密度 ＋ ${name}スコア`:result.source!=="bigwig"?`${result.motif} · モチーフの分布`:`${name} · 信号の分布`;
    $("signal-label").textContent=overlay?"モチーフ密度（左軸）":result.source!=="bigwig"?"モチーフ密度":`${name}（${result.mode==="raw"?"元の値":"平均"}）`;
    $("bw-legend-label").textContent=`${name}（右軸）`;
  }

  function notice(message, type="") {
    $("notice").textContent = message;
    $("notice").className = "notice" + (type ? " " + type : "");
    $("notice").hidden = !message;
  }
  function validInteger(id, min=1) {
    const raw = $(id).value.trim(), n = Number(raw);
    if (!raw || !Number.isSafeInteger(n) || n < min) throw new Error(`${$(id).labels?.[0]?.textContent || id}は${min}以上の整数にしてください。`);
    return n;
  }
  function region() {
    const chrom = $("chrom").value.trim();
    if (!chrom || /\s/.test(chrom)) throw new Error("染色体名を空白を含めずに入力してください。");
    const start = validInteger("start")-1, end = validInteger("end");
    if (end <= start) throw new Error("終了位置は開始位置以上にしてください。");
    return {chrom, start, end};
  }
  async function api(path, data) {
    if (!state.server) throw new Error("この機能にはローカルサーバーが必要です。READMEにある手順で起動し、http://localhost:8765 を開いてください。");
    const response = await fetch(path, {method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify(data)});
    const result = await response.json();
    if (!response.ok) throw new Error(result.error || "サーバーで処理に失敗しました。");
    return result;
  }
  function upload(path, file) {
    if (!state.server) return Promise.reject(new Error("ファイル読み込みにはローカルサーバーが必要です。READMEの起動手順を確認してください。"));
    return new Promise((resolve,reject) => {
      const xhr = new XMLHttpRequest();
      xhr.open("POST",path);
      xhr.setRequestHeader("Content-Type","application/octet-stream");
      xhr.setRequestHeader("X-File-Gzip",/\.gz$/i.test(file.name)?"1":"0");
      xhr.upload.onprogress = event => notice(event.lengthComputable ? `読み込み中: ${file.name} · ${Math.round(event.loaded/event.total*100)}%${event.loaded === event.total ? " · 索引を作成しています…" : ""}` : "読み込み中…");
      xhr.onload = () => {try { const data=JSON.parse(xhr.responseText); if(xhr.status!==200) throw new Error(data.error || "読み込みに失敗しました。"); resolve(data); }catch(e){reject(e);} };
      xhr.onerror = () => reject(new Error("ローカルサーバーへの接続が切れました。サーバーを確認し、読み込み直してください。"));
      xhr.send(file);
    });
  }
  function useBigWig(data,name) {
    dirty();state.bw={...data,name};
    $("bw-info").textContent=`${name} · ${Object.keys(data.chroms).length}染色体`;
    if(!$("sequence").value.trim()&&!state.reference){const [chrom,length]=Object.entries(data.chroms)[0];$("chrom").value=chrom;$("start").value=1;$("end").value=Math.min(length,10000);}
    sync();
  }
  async function importChipAtlas() {
    const accession=$("chip-atlas-id").value.trim().toUpperCase();
    $("chip-atlas-id").value=accession;
    if(!/^[SED]RX\d{6,12}$/.test(accession))throw new Error("SRX・ERX・DRXで始まる実験IDを入力してください（例: SRX000410）。");
    const genome=$("genome-version").value;
    if(genome==="custom")throw new Error("ChIP-Atlasから取得する場合は、対応する参照ゲノム版を選択してください。");
    $("chip-atlas-status").textContent=`${accession}の${genome}版bigWigを取得しています…`;
    notice("ChIP-AtlasからメタデータとbigWigを取得しています。ファイルサイズによって時間がかかります…");
    const data=await api("/api/chip-atlas/import",{accession,genome});
    const details=[data.antigen,data.cell_type].filter(value=>value&&value!=="NA"&&value!=="Unclassified"&&value!=="No description");
    const name=[data.accession,...details].join(" · ");
    useBigWig(data,name);
    $("bw-label").value=data.antigen&&data.antigen!=="NA"?`${data.accession} ${data.antigen}`:data.accession;
    $("chip-atlas-status").textContent=`取得済み · ${data.accession} · ${data.genome} · ${(data.bytes/1024**2).toFixed(1)} MB${data.title?` · ${data.title}`:""}`;
    notice(`ChIP-Atlasから${data.accession}のbigWigを取得しました。解析方法と領域を選んで解析してください。`);
  }
  async function run(task) {
    if(state.busy) return;
    state.busy=true; document.body.classList.add("busy");
    const controls=[...document.querySelectorAll("button,input,textarea,select")].map(e=>[e,e.disabled]);
    controls.forEach(([e])=>e.disabled=true);
    try {await new Promise(r=>setTimeout(r,20)); await task();}
    catch(e) {console.error(e); notice(e.message || String(e),"error");}
    finally {controls.forEach(([e,d])=>e.disabled=d);state.busy=false;document.body.classList.remove("busy");buttons();}
  }
  function buttons() {
    const valid=!!state.result && !state.stale && !state.busy;
    ["save-csv","save-bedgraph","save-bw-bedgraph","save-bed"].forEach(id=>$(id).disabled=!valid);
    $("save-figure").disabled=!state.svg || state.stale || state.busy;
    $("save-guides").disabled=!state.guideResult?.hits.length || state.busy || !state.seq;
    $("table-prev").disabled=!valid || state.page===0;
    $("table-next").disabled=!valid || (state.page+1)*50>=tableData().rows.length;
    $("seq-prev").disabled=state.busy || state.seqPage===0;
    $("seq-next").disabled=state.busy || !state.seq || (state.seqPage+1)*600>=state.seq.sequence.length;
    $("end").readOnly=state.refmode==="paste" && (state.source!=="bigwig" || !!$("sequence").value.trim());
  }
  function dirty(reference=false) {
    if(state.result){state.stale=true;$("chart-subtitle").textContent="条件が変更されています。「解析する」で結果を更新してください。";$("chart-subtitle").classList.add("stale-tag");}
    if(reference){state.seq=null;state.scan=null;state.scanMotif=null;state.guideResult=null;state.svg=null;$("sequence-panel").hidden=true;$("guide-results").hidden=true;$("metric-four").textContent="—";$("guide-status").textContent="参照の変更後は「決定・マーク」で再検索してください。";if(!state.result){$("chart").innerHTML='<div class="empty-state"><h3>参照配列が変更されています</h3><p>「解析する」または「決定・マーク」で表示を更新してください。</p></div>';$("chart-legend").hidden=true;$("chart-range").textContent="参照配列の更新を待機中";}}
    buttons();
  }
  function sync() {
    $("bw-box").hidden=state.source==="sequence";
    $("overlay-hint").hidden=state.source!=="overlay";
    $("ref-optional").hidden=state.source!=="bigwig";
    $("motif-settings").hidden=state.source==="bigwig";
    $("paste-box").hidden=state.refmode!=="paste";
    $("ref-file-box").hidden=state.refmode!=="file"||$("genome-version").value!=="custom";
    $("genome-upload").hidden=$("genome-version").value!=="custom";
    $("location-box").hidden=state.refmode!=="paste";
    $("sliding-settings").hidden=state.mode!=="sliding";
    $("bin-settings").hidden=state.mode!=="bin";
    $("bed-settings").hidden=state.mode!=="bed";
    $("bw-independent-box").hidden=state.source==="sequence";
    const separate=$("bw-independent").checked;
    $("primary-analysis-settings").hidden=state.source==="bigwig"&&separate;
    $("primary-method-label").textContent=state.source==="overlay"&&separate?"モチーフの解析方法":"解析方法";
    $("bw-analysis-settings").hidden=!separate;
    $("bw-window-settings").hidden=$("bw-mode").value!=="sliding";
    $("bw-bin-settings").hidden=$("bw-mode").value!=="bin";
    $("bw-bed-settings").hidden=$("bw-mode").value!=="bed";
    $("coordinate-hint").textContent=state.refmode==="file" ? "参照ゲノムの染色体名と表示領域を指定。BED集計は全参照を対象にします。" : state.source==="bigwig" && !$("sequence").value.trim() ? "bigWigの集計・表示領域を指定。gRNAを表示する場合はこの領域のDNA配列も入力します。" : "貼り付け配列の開始座標を指定します。終了位置は配列長から自動設定。";
    for(const key of ["source","refmode","mode"]){document.querySelectorAll(`[data-${key}]`).forEach(b=>{b.classList.toggle("active",b.dataset[key]===state[key]);b.setAttribute("aria-pressed",String(b.dataset[key]===state[key]));});}
    const f=Number($("flank").value);
    $("window-description").textContent=Number.isSafeInteger(f)&&f>=1?`${fmt(f)} + 中心1 + ${fmt(f)} = ${fmt(2*f+1)} bp / window`:"片側の幅を1以上の整数で指定してください。";
    updateChroms(); buttons();
  }
  function updateChroms() {
    const chroms = new Set([...Object.keys(state.reference?.chroms||{}),...Object.keys(state.bw?.chroms||{})]);
    $("chromosomes").replaceChildren(...[...chroms].map(ch=>{const option=document.createElement("option");option.value=ch;return option;}));
  }
  function pasteLength() {
    try { const seq=C.parseSequence($("sequence").value);$("sequence-length").textContent=`${fmt(seq.length)} bp`;if(state.refmode==="paste"){const start=Number($("start").value);const mapped=state.location;const aligned=mapped&&mapped.query===seq&&mapped.chrom===$("chrom").value.trim()&&mapped.start===start-1;if(Number.isSafeInteger(start)&&start>=1)$("end").value=aligned?mapped.end:start+seq.length-1;} }
    catch {$("sequence-length").textContent=$("sequence").value.trim()?"入力内容を解析時に検証します":"配列を入力してください";}
    sync();
  }
  async function genomeCatalog() {
    if(!state.server)return [];
    const response=await fetch("/api/genomes"),data=await response.json();
    if(!response.ok)throw new Error(data.error||"参照ゲノムの情報を取得できませんでした。");
    return data.genomes;
  }
  function genomeStatus(entry) {
    const names={available:"初回利用時に取得します",queued:"取得の順番を待っています",downloading:"ダウンロード中",verifying:"データを検証中",indexing:"検索用データを準備中",ready:"利用できます",error:"取得に失敗しました"};
    const progress=entry.status==="downloading"&&entry.downloaded?` · ${fmt(Math.round(entry.downloaded/1024**2))} MB${entry.total?` / ${fmt(Math.round(entry.total/1024**2))} MB`:""}`:"";
    return `${entry.name} · ${names[entry.status]||entry.status}${progress}${entry.error?` · ${entry.error}`:""}`;
  }
  async function ensureGenome() {
    const assembly=$("genome-version").value;
    if(assembly==="custom")return;
    if(state.genome?.assembly===assembly&&state.reference?.assembly===assembly)return;
    const name=$("genome-version").selectedOptions[0].textContent;
    notice(`${name}を準備しています。初回の取得・索引作成には数分かかる場合があります。`);
    let data=await api("/api/genomes/select",{assembly});
    while(data.status!=="ready"){
      $("genome-status").textContent=genomeStatus(data);
      if(data.status==="error")throw new Error(data.error||"参照ゲノムの準備に失敗しました。もう一度実行してください。");
      await new Promise(resolve=>setTimeout(resolve,1500));
      const catalog=await genomeCatalog();
      data=catalog.find(entry=>entry.assembly===assembly);
      if(!data)throw new Error("選択した参照ゲノムが見つかりません。");
      if(data.status==="ready")data=await api("/api/genomes/select",{assembly});
    }
    state.genome={...data,name:data.name};state.reference=state.genome;
    $("genome-status").textContent=genomeStatus(data);
    $("genome-info").textContent=`${data.name} · ${Object.keys(data.chroms).length}配列`;
    $("ref-info").textContent=data.name;
    updateChroms();
  }
  async function getSequence(required=true) {
    if(state.refmode==="file")await ensureGenome();
    const r=region();
    if(state.refmode==="paste"){
      if(!$("sequence").value.trim()&&!required)return null;
      let sequence=C.parseSequence($("sequence").value);
      const mapped=state.location;
      if(mapped&&mapped.query===sequence&&mapped.chrom===r.chrom&&mapped.start===r.start){
        if(mapped.referenceSequence)sequence=mapped.referenceSequence;
        else if(mapped.strand==="-")sequence=C.reverseComplement(sequence);
      }
      if(!Number.isSafeInteger(r.start+sequence.length))throw new Error("配列の終端座標が大きすぎます。");
      $("end").value=r.start+sequence.length;
      return {...r,end:r.start+sequence.length,sequence};
    }
    if(!state.reference){if(!required)return null;throw new Error("参照FASTAを選択して読み込んでください。");}
    const data=await api("/api/reference/region",{id:state.reference.id,...r});
    return {...r,sequence:data.sequence};
  }
  async function analyze() {
    if(state.refmode==="file"&&state.source!=="bigwig")await ensureGenome();
    notice("解析しています…");
    if(state.source==="overlay"&&!state.bw)throw new Error("重ねて表示するbigWigファイルを読み込んでください。");
    const motif=state.source!=="bigwig"?C.motifPattern($("motif").value):null;
    const separate=state.source!=="sequence"&&$("bw-independent").checked;
    const usePrimary=state.source!=="bigwig"||!separate;
    const flank=usePrimary&&state.mode==="sliding"?validInteger("flank"):100;
    const step=usePrimary&&state.mode==="sliding"?validInteger("step"):1;
    const bin_size=usePrimary&&state.mode==="bin"?validInteger("bin-size"):null;
    const bed=$("bed-text").value;
    if(usePrimary&&state.mode==="bed")C.parseBed(bed);
    const bwConfig=separate?{mode:$("bw-mode").value,flank:$("bw-mode").value==="sliding"?validInteger("bw-flank"):null,step:$("bw-mode").value==="sliding"?validInteger("bw-step"):null,bed:$("bw-bed-text").value,bin_size:$("bw-mode").value==="bin"?validInteger("bw-bin-size"):null}:{mode:state.mode,flank,step,bed,bin_size};
    if(separate&&bwConfig.mode==="bed")C.parseBed(bwConfig.bed);
    if(!["sliding","bin","bed","raw"].includes(bwConfig.mode))throw new Error("bigWigの解析方法を選択してください。");
    const r=region(), both=$("motif-both").checked;
    let seq=null,scan=null,rows,bwRows=null;
    // File BED can query the whole reference without pulling chromosomes into the browser.
    if(state.source!=="bigwig"&&state.mode==="bed"&&state.refmode==="file"){
      if(!state.reference)throw new Error("参照FASTAを読み込んでください。");
      rows=(await api("/api/reference/bed",{id:state.reference.id,bed,motif,both})).rows;
      if(r.end-r.start<=C.LIMITS.sequence)seq=await getSequence();
      if(seq)scan=C.scanMotif(seq.sequence,motif,both);
    }else if(state.source!=="bigwig"){
      seq=await getSequence();
      const result=C.analyzeSequence(seq.sequence,{chrom:seq.chrom,start:seq.start,mode:state.mode,bed,motif,both,flank,step,binSize:bin_size});
      rows=result.rows;scan=result.scan;
    }else{
      if(!state.bw)throw new Error("解析するbigWigファイルを読み込んでください。");
      // No DNA reference is necessary to summarize a bigWig.
      if(state.refmode==="paste"&&$("sequence").value.trim())seq=await getSequence();
      else if(state.refmode==="file"&&state.reference&&r.end-r.start<=C.LIMITS.sequence)seq=await getSequence();
      rows=(await api("/api/bigwig/analyze",{id:state.bw.id,...bwConfig,...(seq||r),sequence:undefined})).rows;
    }
    if(state.source==="overlay"){
      const signal=(await api("/api/bigwig/analyze",{id:state.bw.id,...bwConfig,...(seq||r),sequence:undefined})).rows;
      bwRows=signal.map(row=>({...row,bw_mean:row.mean,bw_minimum:row.minimum,bw_maximum:row.maximum,bw_coverage:row.coverage,bw_status:row.status}));
      if(!separate){
        if(signal.length!==rows.length)throw new Error("モチーフとbigWigの解析区間が一致していません。");
        rows=rows.map((row,i)=>{
        const bw=signal[i];
        if(row.chrom!==bw.chrom||row.start!==bw.start||row.end!==bw.end)throw new Error("モチーフとbigWigの座標が一致していません。");
        return {...row,bw_mean:bw.mean,bw_minimum:bw.minimum,bw_maximum:bw.maximum,bw_coverage:bw.coverage,bw_status:bw.status};
        });
      }
    }
    state.seq=seq;state.scan=scan;state.scanMotif=motif;state.guideResult=null;state.seqPage=0;
    const mapping=state.refmode==="paste"&&state.location&&state.location.query===$("sequence").value.split(/\r?\n/).filter(line=>!line.trim().startsWith(">")).join("").replace(/\s/g,"").toUpperCase()&&state.location.chrom===(seq||r).chrom&&state.location.start===(seq||r).start?state.location:null;
    state.result={assembly:$("genome-version").value==="custom"?null:$("genome-version").value,rows,bwRows,bwConfig,bin_size:state.source==="bigwig"?bwConfig.bin_size:bin_size,independent:state.source==="overlay"&&separate,source:state.source,mode:state.source==="bigwig"?bwConfig.mode:state.mode,motif,both,flank:state.source==="bigwig"?bwConfig.flank:flank,step:state.source==="bigwig"?bwConfig.step:step,region:seq||r,bwName:state.source==="sequence"?null:state.bw.name,mapping};
    state.stale=false;state.page=0;
    state.tableTrack="motif";
    const available=[...new Set([...rows,...(bwRows||[])].map(row=>row.chrom))];
    state.chartChrom=available.includes(r.chrom)?r.chrom:available[0];
    $("chart-subtitle").classList.remove("stale-tag");
    $("guide-status").textContent="gRNAを入力し、「決定・マーク」で検索してください。";
    render();
    const allRows=state.result.independent?[...rows,...bwRows]:rows;
    const invalid=allRows.filter(row=>row.status!=="ok"||(row.bw_status&&row.bw_status!=="ok")).length;
    const emptySignal=state.source!=="sequence"&&(bwRows||rows).every(row=>row.status!=="ok");
    notice(emptySignal?"bigWigの有効値が指定領域にありません。染色体名・ゲノム版・開始座標を確認してください。座標が不明な場合は「配列を自動検索」を使ってください。":`解析完了 · ${fmt(allRows.length)}区間${state.result.independent?"（両トラック合計）":""}${invalid?` / ${fmt(invalid)}区間は欠測または参照範囲外（表の状態を確認）`:""}${state.mode==="bed"&&state.refmode==="file"&&!seq?"。表示領域が2 Mbを超えるため配列ビューは省略しています。":""}`,invalid?"warn":"");
  }
  async function applyLocation(match,query){
    if(!match||query!==C.parseSequence($("sequence").value))throw new Error("検索配列が変わっています。もう一度自動検索してください。");
    if(state.bw&&(!(match.chrom in state.bw.chroms)||match.end>state.bw.chroms[match.chrom]))throw new Error("一致した座標がbigWigの染色体・範囲と合いません。検索用参照とbigWigのゲノム版を確認してください。");
    if(state.bw&&!match.record.includes(":")&&state.genome.chroms[match.record]!==state.bw.chroms[match.chrom])throw new Error("検索用参照とbigWigの染色体長が違います。同じゲノム版のデータを使用してください。");
    let referenceSequence=null;
    if(match.edits)referenceSequence=(await api("/api/reference/region",{id:state.genome.id,chrom:match.record,start:match.local_start,end:match.local_end})).sequence;
    dirty(true);state.location={...match,query,referenceSequence};
    $("chrom").value=match.chrom;$("start").value=match.start+1;$("end").value=match.end;
    $("location-status").textContent=`一致: ${match.chrom}:${fmt(match.start+1)}–${fmt(match.end)} (${match.strand}鎖)${match.edits?" · 入力と1塩基差あり。参照ゲノム配列を使用して解析":match.strand==="-"?" · 解析はゲノム＋鎖の向きに変換":""}`;
    if(state.bw)state.source="overlay";
    pasteLength();await analyze();
    if($("guides").value.trim())await markGuides();
  }
  async function locate(){
    const query=C.parseSequence($("sequence").value);
    if(query.length<20)throw new Error("位置検索には20 bp以上の配列を入力してください。");
    await ensureGenome();
    if(!state.genome)throw new Error("「自分のFASTAを使用」では検索用参照FASTAを読み込んでください。");
    notice("参照ゲノム内を検索しています…大きなゲノムでは時間がかかります。");
    const data=await api("/api/reference/search",{id:state.genome.id,query,max_edits:1});
    state.locationCandidates=data.matches.map(match=>({...match,query}));
    $("location-results").hidden=true;
    if(!data.matches.length)throw new Error("完全一致または1塩基差以内の位置が見つかりません。参照ゲノムの版・配列の違い・検索対象範囲を確認してください。");
    if(data.truncated)throw new Error("一致が50か所を超えています。より長く特徴的な配列を使用してください。座標は自動設定していません。");
    if(data.matches.length===1&&!data.matches[0].edits)return applyLocation(data.matches[0],query);
    const differences={query_deletion:"入力に1塩基欠失",query_insertion:"入力に1塩基挿入",substitution:"1塩基置換",exact:"完全一致"};
    $("location-match").replaceChildren(...data.matches.map((match,index)=>{const option=document.createElement("option");option.value=index;option.textContent=`${match.chrom}:${fmt(match.start+1)}–${fmt(match.end)} (${match.strand}) · ${differences[match.difference]||"完全一致"}`;return option;}));
    $("location-match").value="0";$("location-results").hidden=false;
    const hasDifference=data.matches.some(match=>match.edits);
    $("location-status").textContent=`${data.matches.length}候補。${hasDifference?"1塩基差があります。選択後は参照ゲノム配列で解析します。入力配列は書き換えません。":"使用する座標を選んでください。"}`;
    $("apply-location").textContent=hasDifference?"この候補の参照配列で解析":"この座標で解析";
    notice(hasDifference?"入力配列と1塩基違う候補が見つかりました。差と座標を確認して選択してください。選択後は参照ゲノム配列でCG密度・gRNAを解析します。":"複数の位置に一致しました。候補を選んで「この座標で解析」を押してください。","warn");
  }
  async function markGuides() {
    const seq=await getSequence();
    // Reference edits invalidate analysis through dirty(); don't combine mismatched coordinates.
    const result=C.findGuides(seq.sequence,$("guides").value,{chrom:seq.chrom,start:seq.start,both:$("guide-both").checked});
    state.seq=seq;state.guideResult=result;state.seqPage=0;
    if(state.source!=="bigwig"){state.scanMotif=C.motifPattern($("motif").value);state.scan=C.scanMotif(seq.sequence,state.scanMotif,$("motif-both").checked);}
    else {state.scan=null;state.scanMotif=null;}
    renderGuides();renderSequence();renderChart();buttons();
    $("metric-four").textContent=fmt(result.hits.length);
    $("guide-status").textContent=`${result.guides.length}本 / ${fmt(result.hits.length)}一致 · ${seq.chrom}:${fmt(seq.start+1)}–${fmt(seq.end)} 内を検索`;
    notice(result.unmatched.length?`gRNAをマークしました。一致なし: ${result.unmatched.join("、")}（表示領域内を完全一致で検索）`:`gRNAをマークしました · ${fmt(result.hits.length)}一致`,result.unmatched.length?"warn":"");
  }
  function render() {
    const result=state.result,rows=result.rows,valid=rows.filter(r=>r.status==="ok");
    $("metric-one-label").textContent=result.source!=="bigwig"?"参照配列長":"bigWigデータ";
    $("metric-one").textContent=result.source!=="bigwig"?fmt(state.refmode==="file"?Object.values(state.reference.chroms).reduce((a,b)=>a+b,0):(state.seq?.sequence.length||0)):fmt(Object.keys(state.bw.chroms).length);
    $("metric-one-unit").textContent=result.source!=="bigwig"?"bp":"chromosomes";
    $("metric-two-label").textContent=result.source!=="bigwig"?"表示領域のモチーフ一致":result.mode==="raw"?"区間値の単純平均":"有効区間の平均信号";
    $("metric-two").textContent=result.source!=="bigwig"?(state.scan?fmt(state.scan.prefix.at(-1)):"—"):(valid.length?number(valid.reduce((n,r)=>n+r.mean,0)/valid.length):"—");
    $("metric-two-unit").textContent=result.source!=="bigwig"?"sites":"mean of means";
    const overlay=result.source==="overlay";
    $("metric-three-label").textContent=overlay?"有効区間：モチーフ / bigWig":"有効な解析区間";
    $("metric-three").textContent=overlay?`${fmt(valid.length)} / ${fmt((result.bwRows||rows).filter(r=>(r.bw_status||r.status)==="ok").length)}`:fmt(valid.length);
    $("metric-three-unit").textContent="regions";$("metric-four").textContent="—";
    refreshFigureLabels(result);
    $("chart-subtitle").textContent=result.independent?`モチーフ: ${configText(result)} ／ ${bwDisplayName()}: ${configText(result.bwConfig)}`:configText(result);
    if(result.assembly)$("chart-subtitle").textContent+=` · 参照ゲノム: ${result.assembly}`;
    if(result.mapping?.edits)$("chart-subtitle").textContent+=" · 入力と1塩基差あり／参照ゲノム配列で解析";
    $("bw-legend").hidden=!overlay;
    $("save-bw-bedgraph").hidden=!overlay;
    $("save-bedgraph").textContent=overlay?"モチーフ bedGraph ↓":"bedGraph ↓";
    $("chart-note").textContent=result.independent?"左右の縦軸は別スケール · 各トラックの解析方法は上部に表示":overlay?"左右の縦軸は別スケール · bigWigは同じ区間の平均値":result.source!=="bigwig"?"密度 = モチーフ開始位置の数 / 区間長":result.mode==="raw"?"元のbigWig値 · 追加の平滑化なし":"平均信号 = データがある塩基の平均（欠測は除外）";
    $("result-track-tabs").hidden=!result.independent;
    $("table-summary").textContent=`全${fmt(rows.length)}件 · 座標は1-based、両端を含む · CSVには欠測を含む全行を保存`;
    renderChart();renderSequence();renderTable();renderGuides();buttons();
  }
  function renderTable() {
    if(!state.result)return;
    const {rows,source,mode}=tableData();
    $("table-motif").setAttribute("aria-pressed",String(state.tableTrack!=="bigwig"));
    $("table-bigwig").setAttribute("aria-pressed",String(state.tableTrack==="bigwig"));
    $("table-motif").classList.toggle("active",state.tableTrack!=="bigwig");
    $("table-bigwig").classList.toggle("active",state.tableTrack==="bigwig");
    $("save-bed").textContent=state.result.independent?"表示トラックのBED ↓":"BED ↓";
    $("table-summary").textContent=state.result.independent?`${source==="bigwig"?bwDisplayName():"モチーフ"}: ${fmt(rows.length)}件 · CSVは両トラックを別行に保存`:`全${fmt(rows.length)}件 · 座標は1-based、両端を含む · CSVには欠測を含む全行を保存`;
    const fields=source==="overlay"?[["count","モチーフ数"],["density","密度"],["bw_mean","bigWig平均"],["bw_minimum","bigWig最小"],["bw_maximum","bigWig最大"],["bw_coverage","bigWig被覆率"]]:source!=="bigwig"?[["count","モチーフ数"],["density","密度"]]:mode==="raw"?[["mean","元のbigWig値"]]:[["mean","平均"],["minimum","最小"],["maximum","最大"],["coverage","被覆率"]];
    const selected=rows.slice(state.page*50,(state.page+1)*50);
    const status=(s)=>`<span class="status ${s==="ok"?"":"missing"}">${esc(statuses[s]||s)}</span>`;
    $("table-wrap").innerHTML=`<table><thead><tr><th>領域名</th><th>染色体</th><th>開始</th><th>終了</th>${fields.map(([,l])=>`<th>${l}</th>`).join("")}<th>${source==="overlay"?"モチーフ状態":"状態"}</th>${source==="overlay"?"<th>bigWig状態</th>":""}</tr></thead><tbody>${selected.map(r=>`<tr><td>${esc(r.name)}</td><td>${esc(r.chrom)}</td><td class="numeric">${fmt(r.start+1)}</td><td class="numeric">${fmt(r.end)}</td>${fields.map(([k])=>`<td class="numeric">${k==="count"?(r[k]==null?"—":fmt(r[k])):k.endsWith("coverage")?(r[k]==null?"—":(r[k]*100).toFixed(1)+"%"):esc(number(r[k]))}</td>`).join("")}<td>${status(r.status)}</td>${source==="overlay"?`<td>${status(r.bw_status)}</td>`:""}</tr>`).join("")}</tbody></table>`;
    $("table-page-info").textContent=`${fmt(state.page*50+1)}–${fmt(Math.min((state.page+1)*50,rows.length))} / ${fmt(rows.length)}件`;
    buttons();
  }
  function renderGuides() {
    const result=state.guideResult;
    $("guide-results").hidden=!result;
    if(!result)return;
    const chrom=state.seq.chrom;
    $("guide-summary").textContent=`${result.guides.length}本 / ${fmt(result.hits.length)}一致${result.unmatched.length?` · 一致なし: ${result.unmatched.join("、")}`:""} · 表は先頭200一致`;
    $("guide-table").innerHTML=`<table><thead><tr><th>名前</th><th>染色体</th><th>開始</th><th>終了</th><th>鎖</th><th>配列（入力5′→3′）</th></tr></thead><tbody>${result.hits.slice(0,200).map(h=>`<tr><td>${esc(h.name)}</td><td>${esc(chrom)}</td><td class="numeric">${fmt(h.start+1)}</td><td class="numeric">${fmt(h.end)}</td><td>${h.strand}</td><td>${esc(h.sequence)}</td></tr>`).join("")||'<tr><td colspan="6">表示領域内に一致はありません。</td></tr>'}</tbody></table>`;
  }
  function renderSequence() {
    const seq=state.seq;
    $("sequence-panel").hidden=!seq;
    if(!seq)return;
    const start=state.seqPage*600,end=Math.min(start+600,seq.sequence.length);
    const motifMask=new Uint8Array(end-start),guideMask=new Uint8Array(end-start);
    if(state.scan){const length=state.scanMotif?.length||1;for(let i=Math.max(0,start-length+1);i<end;i++){if(state.scan.flags[i])for(let j=Math.max(i,start);j<Math.min(i+length,end);j++)motifMask[j-start]=1;}}
    for(const hit of state.guideResult?.hits||[]){for(let j=Math.max(start,hit.start-seq.start);j<Math.min(end,hit.end-seq.start);j++)guideMask[j-start]=1;}
    let html="";
    for(let i=start;i<end;i+=60){html+=`<div class="seq-line"><span class="seq-position">${fmt(seq.start+i+1)}</span><span class="bases">`;for(let j=i;j<Math.min(i+60,end);j++){const b=seq.sequence[j];html+=`<span class="base${motifMask[j-start]?" motif":""}${guideMask[j-start]?" guide":""}${/[^ACGT]/.test(b)?" unknown":""}" title="${esc(seq.chrom)}:${seq.start+j+1}">${b}</span>${(j-i+1)%10===0?" ":""}`;}html+="</span></div>";}
    $("sequence-view").innerHTML=html;
    $("sequence-caption").textContent=`${seq.chrom}:${fmt(seq.start+start+1)}–${fmt(seq.start+end)} · マークは完全一致の全長`;
    buttons();
  }
  function renderChart() {
    const result=state.stale?null:state.result,seq=state.seq;
    if(!result&&!seq)return;
    const chrom=result?state.chartChrom:seq?.chrom;
    const rows=(result?.rows||[]).filter(r=>r.chrom===chrom),source=result?.source||'sequence',overlay=source==='overlay';
    const signalRows=(result?.bwRows||result?.rows||[]).filter(r=>r.chrom===chrom);
    let left=Infinity,right=0;
    const domainRows=overlay?[...rows,...signalRows]:rows;
    for(const row of domainRows){left=Math.min(left,row.start);right=Math.max(right,row.end);}
    if(!domainRows.length){left=seq.start;right=seq.end;}
    if(right<=left)right=left+1;
    const showConfig=result?.independent||result?.mode==='bin';
    const W=1040,L=75,R=overlay?90:25,T=showConfig?58:40,B=290;
    const x=p=>L+(p-left)/(right-left)*(W-L-R);
    function makeTrack(field,statusField,color,label,dashed=false,trackRows=rows,mode=result?.mode){
      const valid=trackRows.filter(r=>r[statusField]==='ok'&&Number.isFinite(r[field]));
      let min=0,max=0;for(const row of valid){min=Math.min(min,row[field]);max=Math.max(max,row[field]);}
      if(max===min)max=min+1;const pad=(max-min)*.12;max+=pad;if(min<0)min-=pad;
      return {field,statusField,color,label,dashed,valid,rows:trackRows,mode,min,max,y:v=>B-(v-min)/(max-min)*(B-T)};
    }
    const bwLabel=bwDisplayName();
    const primary=makeTrack(source==='bigwig'?'mean':'density','status','#147d77',source==='bigwig'?bwLabel+(result?.mode==='raw'?' value':' mean'):'Motif count / bp');
    const secondary=overlay?makeTrack('bw_mean','bw_status','#8b5db7',bwLabel+(result?.bwConfig?.mode==='raw'?' value':' mean'),true,signalRows,result?.bwConfig?.mode||result?.mode):null;
    const tracks=secondary?[primary,secondary]:[primary];
    const guideHits=(state.guideResult?.hits||[]).filter(h=>seq?.chrom===chrom&&h.end>left&&h.start<right);
    const palette=['#cb894b','#9b63bc','#d84c60','#118c75','#946536','#3888b4','#b25887','#718c2f'];
    const guideKey=guide=>guide.name+'\u0000'+guide.sequence;
    const colors=new Map((state.guideResult?.guides||[]).map((guide,index)=>[guideKey(guide),palette[index]||`hsl(${Math.round(index*137.508)%360},60%,40%)`]));
    const laneEnds=[],guideMarkers=[];
    let hiddenLabels=0;
    for(const hit of guideHits.slice(0,500)){
      // Midpoint of the first and last matched bases; this is not a cut-site estimate.
      const midpoint=(hit.start+hit.end-1)/2,mx=x(midpoint);
      if(midpoint<left||midpoint>=right)continue;
      const name=Array.from(hit.name).length>30?Array.from(hit.name).slice(0,29).join('')+'…':hit.name;
      const range=`${fmt(hit.start+1)}–${fmt(hit.end)} (${hit.strand})`;
      const nameWidth=Array.from(name).reduce((n,c)=>n+(/[^\x00-\x7f]/.test(c)?11:6.5),0);
      const width=Math.min(W-L-R,Math.max(90,nameWidth,range.length*6)+14);
      const labelX=Math.max(L+width/2,Math.min(W-R-width/2,mx));
      let lane=laneEnds.findIndex(end=>labelX-width/2>end+10);
      if(lane<0&&laneEnds.length<8){lane=laneEnds.length;laneEnds.push(0);}
      if(lane>=0)laneEnds[lane]=labelX+width/2;else hiddenLabels++;
      guideMarkers.push({hit,mx,name,range,labelX,lane,color:colors.get(guideKey(hit))||palette[0]});
    }
    const guideTop=B+58;
    const H=Math.max(430,guideTop+Math.max(0,laneEnds.length-1)*46+82);
    const chroms=[...new Set(result?[...result.rows,...(result.bwRows||[])].map(r=>r.chrom):[chrom])];
    let svg=`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(chrom)}の${esc(overlay?'モチーフ密度と'+bwLabel+'スコアの重ね合わせ':source==='bigwig'?bwLabel+'平均信号':'モチーフ密度')}とgRNA位置"><desc>${overlay?'緑実線は左軸のモチーフ密度、紫破線は右軸の'+esc(bwLabel)+'平均スコア。左右の軸は別スケール。':''}欠測は0にせず線を切ります。</desc><rect width="${W}" height="${H}" fill="white"/><defs><clipPath id="plot-clip"><rect x="${L}" y="${T}" width="${W-L-R}" height="${B-T}"/></clipPath></defs>`;
    if(result?.assembly)svg+=`<text x="${L}" y="${H-43}" fill="#718488" font-size="10" font-family="sans-serif">Reference genome: ${esc(result.assembly)}</text>`;
    svg+=`<line x1="${L}" x2="${L+22}" y1="15" y2="15" stroke="${primary.color}" stroke-width="2"/><text x="${L+29}" y="19" fill="${primary.color}" font-size="11" font-family="sans-serif">${esc(overlay?'Motif density (left)':primary.label)}</text>`;
    if(secondary)svg+=`<line x1="${L+230}" x2="${L+252}" y1="15" y2="15" stroke="${secondary.color}" stroke-width="2" stroke-dasharray="6 3"/><text x="${L+259}" y="19" fill="${secondary.color}" font-size="11" font-family="sans-serif">${esc(bwLabel)} score (right)</text>`;
    if(result?.independent){const brief=c=>c.mode==='raw'?'stored values (no smoothing)':c.mode==='bed'?'BED':c.mode==='bin'?`${c.bin_size} bp bins (last bin clipped)`:`${2*c.flank+1} bp window / ${c.step} bp step`;svg+=`<text x="${L}" y="34" fill="#718488" font-size="10" font-family="sans-serif">Motif: ${brief(result)} | ${esc(bwLabel)}: ${brief(result.bwConfig)}</text>`;}
    if(result?.mode==='bin'&&!result.independent)svg+=`<text x="${L}" y="34" fill="#718488" font-size="10" font-family="sans-serif">${result.bin_size} bp bins, anchored at region start; last bin clipped</text>`;
    if(result?.mapping?.edits)svg+=`<text x="${L}" y="${showConfig?50:34}" fill="#966631" font-size="10" font-family="sans-serif">Reference-genome sequence used (query differs by 1 base)</text>`;
    for(let i=0;i<=4;i++){
      const v=primary.min+(primary.max-primary.min)*i/4,yy=primary.y(v);
      svg+=`<line x1="${L}" y1="${yy}" x2="${W-R}" y2="${yy}" stroke="#e6eeea"/><text x="${L-12}" y="${yy+4}" text-anchor="end" fill="${overlay?primary.color:'#718488'}" font-size="11" font-family="sans-serif">${esc(number(v))}</text>`;
      if(secondary){const v2=secondary.min+(secondary.max-secondary.min)*i/4;svg+=`<text x="${W-R+12}" y="${yy+4}" fill="${secondary.color}" font-size="11" font-family="sans-serif">${esc(number(v2))}</text>`;}
    }
    for(let i=0;i<=5;i++){const p=left+(right-left)*i/5,xx=x(p);svg+=`<text x="${xx}" y="${B+22}" text-anchor="middle" fill="#718488" font-size="11" font-family="sans-serif">${fmt(Math.min(right,Math.floor(p)+1))}</text>`;}
    svg+=`<text x="18" y="${(T+B)/2}" transform="rotate(-90 18 ${(T+B)/2})" text-anchor="middle" font-size="11" fill="${overlay?primary.color:'#718488'}" font-family="sans-serif">${esc(axisLabel(primary.label))}<title>${esc(primary.label)}</title></text>`;
    if(secondary)svg+=`<text x="${W-18}" y="${(T+B)/2}" transform="rotate(90 ${W-18} ${(T+B)/2})" text-anchor="middle" font-size="11" fill="${secondary.color}" font-family="sans-serif">${esc(axisLabel(secondary.label))}<title>${esc(secondary.label)}</title></text>`;
    function drawTrack(track){
      const {field,statusField,color,y,dashed}=track,dash=dashed?' stroke-dasharray="6 3"':'';
      const isValid=r=>r[statusField]==='ok'&&Number.isFinite(r[field]);
      let out=`<g data-track="${field}" data-method="${esc(track.mode||'none')}" clip-path="url(#plot-clip)">`;
      if(track.mode==='raw'){
        let path='',previousEnd=null;
        for(const row of track.rows){if(!isValid(row)){previousEnd=null;continue;}path+=`${previousEnd===row.start?'L':'M'}${x(row.start).toFixed(2)},${y(row[field]).toFixed(2)}H${x(row.end).toFixed(2)} `;previousEnd=row.end;}
        out+=`<path d="${path}" fill="none" stroke="${color}" stroke-width="2"${dash}/>`;
      }else if(track.mode==='bed'||track.mode==='bin'){
        const hop=Math.max(1,Math.ceil(track.valid.length/3000));
        for(let i=0;i<track.valid.length;i+=hop){const r=track.valid[i],xx=x(r.start),yy=y(r[field]);out+=`<rect x="${xx}" y="${Math.min(y(0),yy)}" width="${Math.max(1,x(r.end)-xx)}" height="${Math.max(1,Math.abs(yy-y(0)))}" fill="${color}" opacity=".10"><title>${esc(r.name)} · ${esc(chrom)}:${r.start+1}–${r.end}\n${esc(track.label)}: ${number(r[field])}</title></rect><line x1="${xx}" x2="${x(r.end)}" y1="${yy}" y2="${yy}" stroke="${color}" stroke-width="2"${dash}><title>${esc(track.label)}: ${number(r[field])}</title></line>`;}
      }else{
        const plotted=[],bucket=Math.max(1,Math.ceil(track.rows.length/1200));
        for(let i=0;i<track.rows.length;i+=bucket){let min=null,max=null;const flush=()=>{if(min&&max){if(min.start<=max.start){plotted.push(min);if(max!==min)plotted.push(max);}else{plotted.push(max);plotted.push(min);}}min=max=null;};for(const r of track.rows.slice(i,i+bucket)){if(!isValid(r)){flush();plotted.push(r);continue;}if(!min||r[field]<min[field])min=r;if(!max||r[field]>max[field])max=r;}flush();}
        let path='',pen=false,segmentCount=0,segmentPoint=null,dots='';
        const finishSegment=()=>{if(segmentCount===1)dots+=`<circle cx="${x(segmentPoint.center??(segmentPoint.start+segmentPoint.end)/2)}" cy="${y(segmentPoint[field])}" r="3" fill="${color}"/>`;segmentCount=0;segmentPoint=null;};
        for(const r of plotted){if(!isValid(r)){finishSegment();pen=false;continue;}path+=`${pen?'L':'M'}${x(r.center??(r.start+r.end)/2).toFixed(2)},${y(r[field]).toFixed(2)} `;pen=true;segmentCount++;segmentPoint=r;}
        finishSegment();
        out+=`<path d="${path}" fill="none" stroke="${color}" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"${dash}/>`;
        out+=dots;
      }
      return out+'</g>';
    }
    svg+='<g data-guide-overlay="true" clip-path="url(#plot-clip)">';
    for(const marker of guideMarkers){const {hit,mx,color}=marker,leftEdge=x(hit.start-.5),rightEdge=x(hit.end-.5);svg+=`<rect x="${leftEdge}" y="${T}" width="${Math.max(2,rightEdge-leftEdge)}" height="${B-T}" fill="${color}" opacity=".09"/><line data-guide-midpoint="${(hit.start+hit.end-1)/2}" x1="${mx}" x2="${mx}" y1="${T}" y2="${B}" stroke="${color}" stroke-width="1.3" stroke-dasharray="4 4" opacity=".85"><title>${esc(hit.name)} · ${esc(chrom)}:${hit.start+1}–${hit.end} (${hit.strand}) · sequence midpoint</title></line>`;}
    svg+='</g>';
    for(const track of tracks)svg+=drawTrack(track);
    if(!tracks.some(track=>track.valid.length))svg+=`<text x="${W/2}" y="${(T+B)/2}" text-anchor="middle" fill="#8b9d98" font-size="13" font-family="sans-serif">この染色体に有効な信号値はありません</text>`;
    else if(secondary&&!secondary.valid.length)svg+=`<text x="${L+8}" y="${T+20}" fill="${secondary.color}" font-size="12" font-family="sans-serif">${esc(bwLabel)}: この領域にデータがありません。染色体名・開始座標を確認してください。</text>`;
    svg+='<g data-guide-labels="true">';
    for(const marker of guideMarkers){const {hit,mx,name,range,labelX,lane,color}=marker;if(lane<0)continue;const yy=guideTop+lane*46;svg+=`<text data-guide-label="${esc(hit.name)}" x="${labelX}" y="${yy}" text-anchor="middle" fill="${color}" font-size="11" font-weight="bold" font-family="sans-serif">${esc(name)}<title>${esc(hit.name)}</title></text><text x="${labelX}" y="${yy+13}" text-anchor="middle" fill="${color}" font-size="10" font-family="sans-serif">${esc(range)}</text><line x1="${x(Math.max(left,hit.start))}" x2="${x(Math.min(right,hit.end-1))}" y1="${yy+27}" y2="${yy+27}" stroke="${color}" stroke-width="5" stroke-linecap="round"/><circle cx="${mx}" cy="${yy+27}" r="3.6" fill="${color}"/>`;}
    svg+='</g>';
    svg+=`<text x="${L}" y="${H-28}" fill="#8b9c98" font-size="10" font-family="sans-serif">gRNA · dashed lines: sequence midpoint (not cut sites)${guideHits.length>500?' · first 500 matches shown':''}${hiddenLabels?' · dense labels omitted; see table':''}${guideHits.length?'':' · no marks'}</text><text x="${W/2}" y="${H-6}" text-anchor="middle" fill="#718488" font-size="11" font-family="sans-serif">${esc(chrom)} · Position (1-based)</text></svg>`;
    state.svg=svg;
    $('chart').innerHTML=(chroms.length>1?`<div class="result-chrom"><label for="chart-chrom">表示する染色体</label><select id="chart-chrom">${chroms.map(ch=>`<option${ch===chrom?' selected':''}>${esc(ch)}</option>`).join('')}</select></div>`:'')+svg;
    $('chart-chrom')?.addEventListener('change',event=>{state.chartChrom=event.target.value;renderChart();});
    $('chart-legend').hidden=false;$('bw-legend').hidden=!overlay;
    $('chart-range').textContent=`${chrom}:${fmt(left+1)}–${fmt(right)} · ${result?.independent?`モチーフ ${fmt(rows.length)} / bigWig ${fmt(signalRows.length)}`:fmt(rows.length)}区間${rows.length>3000?' · 描画のみ間引きあり':''}`;
    const chartSvg=$('chart').querySelector('svg'),tooltip=document.createElement('div');tooltip.className='chart-tooltip';tooltip.hidden=true;$('chart').append(tooltip);
    if(result&&domainRows.length){chartSvg.addEventListener('pointermove',event=>{const box=chartSvg.getBoundingClientRect(),xx=(event.clientX-box.left)/box.width*W,yy=(event.clientY-box.top)/box.height*H;if(xx<L||xx>W-R||yy<T||yy>B){tooltip.hidden=true;return;}const pos=left+(xx-L)/(W-L-R)*(right-left);tooltip.textContent=`${chrom}:${fmt(Math.floor(pos)+1)}\n`+tracks.map(track=>{const list=track.rows;let row=null;if(track.mode==='sliding'&&list.length){let lo=0,hi=list.length-1;while(lo<hi){const mid=Math.floor((lo+hi)/2);if(list[mid].center<pos)lo=mid+1;else hi=mid;}row=list[lo];if(lo>0&&Math.abs(list[lo-1].center-pos)<Math.abs(row.center-pos))row=list[lo-1];}else row=list.find(r=>r.start<=pos&&pos<r.end);return `${track.label}: ${!row?'対象区間なし':row[track.statusField]==='ok'?number(row[track.field]):statuses[row[track.statusField]]}${row?` (${fmt(row.start+1)}–${fmt(row.end)})`:''}`;}).join('\n');tooltip.hidden=false;tooltip.style.left=Math.max(0,Math.min(event.clientX-box.left+12,box.width-310))+'px';tooltip.style.top=Math.max(0,event.clientY-box.top-30)+'px';});chartSvg.addEventListener('pointerleave',()=>tooltip.hidden=true);}
    buttons();
  }
  function download(text, filename, type="text/plain;charset=utf-8") {
    downloadBlob(new Blob([text],{type}),filename);
  }
  function downloadBlob(blob,filename) {
    const url=URL.createObjectURL(blob),a=document.createElement("a");a.href=url;a.download=filename;document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),2000);
  }
  async function saveFigure() {
    if(!state.svg||state.stale)return;
    if($("figure-format").value==="svg"){download(state.svg,"motif_studio.svg","image/svg+xml;charset=utf-8");return;}
    const viewBox=state.svg.match(/viewBox="[\d.-]+ [\d.-]+ ([\d.]+) ([\d.]+)"/);
    if(!viewBox)throw new Error("図の大きさを取得できませんでした。");
    const width=Number(viewBox[1]),height=Number(viewBox[2]),scale=2;
    const source=URL.createObjectURL(new Blob([state.svg],{type:"image/svg+xml;charset=utf-8"}));
    try {
      const image=new Image();
      await new Promise((resolve,reject)=>{image.onload=resolve;image.onerror=()=>reject(new Error("PNG用の図を作成できませんでした。"));image.src=source;});
      const canvas=document.createElement("canvas");canvas.width=Math.round(width*scale);canvas.height=Math.round(height*scale);
      const context=canvas.getContext("2d");if(!context)throw new Error("このブラウザではPNGを作成できません。");
      context.fillStyle="#ffffff";context.fillRect(0,0,canvas.width,canvas.height);context.drawImage(image,0,0,canvas.width,canvas.height);
      const blob=await new Promise(resolve=>canvas.toBlob(resolve,"image/png"));
      if(!blob)throw new Error("PNGの変換に失敗しました。");
      downloadBlob(blob,"motif_studio.png");
    } finally {URL.revokeObjectURL(source);}
  }
  function safeName(name){return String(name).replace(/[^a-zA-Z0-9_.-]/g,"_").slice(0,80)||"result";}
  function exportResults(kind) {
    if(!state.result||state.stale)return;
    const result=state.result,overlay=result.source==="overlay",bwTrack=kind==="bw-bedgraph";
    const useBwRows=!!(bwTrack&&result.bwRows||result.independent&&kind==="bed"&&state.tableTrack==="bigwig");
    const sourceRows=useBwRows?result.bwRows:result.rows;
    const config=useBwRows?result.bwConfig:result;
    const field=useBwRows?"mean":bwTrack?"bw_mean":result.source!=="bigwig"?"density":"mean";
    const base=`${safeName(useBwRows||bwTrack?"bigwig":result.motif||"bigwig")}_${config.mode}`;
    if(kind==="csv"){
      const keys=["assembly",...(result.mode==="bin"?["bin_bp","interval_bp"]:[]),"chrom","start","end","name","center",...(result.source!=="bigwig"?["count","density"]:["mean","minimum","maximum","coverage"]),"status",...(overlay?["bw_mean","bw_minimum","bw_maximum","bw_coverage","bw_status"]:[]),...(result.mapping?["sequence_basis","query_difference","query_strand"]:[])];
      const csvCell=value=>'"'+String(value??"").replace(/"/g,'""')+'"';
      const userCell=value=>/^[=+\-@\t\r]/.test(String(value))?"'"+value:value;
      if(result.independent){
        const columns=["assembly","track","track_label","chrom","start_0based","end_exclusive","center_0based","name","method","flank_bp","step_bp",...([result.mode,result.bwConfig.mode].includes("bin")?["bin_bp","interval_bp"]:[]),"count","density","mean","minimum","maximum","coverage","status","sequence_basis","query_difference","query_strand"];
        const lines=[columns.join(",")];
        for(const [track,list,cfg] of [["motif",result.rows,result],["bigwig",result.bwRows,result.bwConfig]])for(const row of list){const data={...row,assembly:result.assembly,track,track_label:track==="motif"?result.motif:bwDisplayName(),start_0based:row.start,end_exclusive:row.end,center_0based:row.center,method:cfg.mode,flank_bp:cfg.mode==="sliding"?cfg.flank:null,step_bp:cfg.mode==="sliding"?cfg.step:null,bin_bp:cfg.mode==="bin"?cfg.bin_size:null,interval_bp:row.end-row.start,sequence_basis:track==="bigwig"?"stored_bigwig":result.mapping?.edits?"reference_genome":"query_genomic_forward",query_difference:result.mapping?.difference,query_strand:result.mapping?.strand};lines.push(columns.map(key=>csvCell(["name","chrom","track_label"].includes(key)?userCell(data[key]):data[key])).join(","));}
        download("\uFEFF"+lines.join("\r\n"),`${safeName(result.motif)}_and_bigwig_independent.csv`,"text/csv;charset=utf-8");return;
      }
      const lines=[keys.map(k=>k==="start"?"start_0based":k==="end"?"end_exclusive":k==="center"?"center_0based":k).join(",")];
      for(const row of result.rows){const data={...row,assembly:result.assembly,bin_bp:result.bin_size,interval_bp:row.end-row.start,sequence_basis:result.mapping?.edits?"reference_genome":"query_genomic_forward",query_difference:result.mapping?.difference,query_strand:result.mapping?.strand};lines.push(keys.map(k=>csvCell(["name","chrom"].includes(k)?userCell(data[k]):data[k])).join(","));}
      download("\uFEFF"+lines.join("\r\n"),base+".csv","text/csv;charset=utf-8");return;
    }
    const rows=sourceRows.filter(r=>kind==="bed"?r.status==="ok"||(!result.independent&&overlay&&r.bw_status==="ok"):r[useBwRows?"status":bwTrack?"bw_status":"status"]==="ok"&&Number.isFinite(r[field])).slice().sort((a,b)=>a.chrom.localeCompare(b.chrom)||a.start-b.start||a.end-b.end);
    if(!rows.length){notice("このトラックに出力できる有効な区間がありません。CSVで欠測の理由を確認してください。","warn");return;}
    if(kind==="bed"){
      download(rows.map(r=>[r.chrom,r.start,r.end,safeName(r.name),0,"."].join("\t")).join("\n")+"\n",base+".bed");return;
    }
    // Signal windows overlap. A bedGraph track must be non-overlapping: sliding uses center points.
    const output=rows.map(r=>({...r,start:config.mode==="sliding"?r.center:r.start,end:config.mode==="sliding"?r.center+1:r.end}));
    for(let i=1;i<output.length;i++)if(output[i].chrom===output[i-1].chrom&&output[i].start<output[i-1].end){notice("BED領域が重複しているためbedGraphには出力できません。CSVまたはBEDで出力するか、重複しない領域を指定してください。","warn");return;}
    download(output.map(r=>[r.chrom,r.start,r.end,r[field]].join("\t")).join("\n")+"\n",base+".bedGraph");
  }
  function demo() {
    let seed=42;const random=()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return seed/4294967296;};
    const bases="ATGC";let seq="";for(let i=0;i<2400;i++){seq+=random()<(i>650&&i<1600?.22:.025)?"CG":bases[Math.floor(random()*4)];}seq=seq.slice(0,2400);
    state.source="sequence";state.refmode="paste";state.mode="sliding";state.result=null;state.stale=false;state.seq=null;state.scan=null;state.guideResult=null;state.chartChrom=null;
    $("sequence").value=">chrDemo\n"+seq.match(/.{1,80}/g).join("\n");$("chrom").value="chrDemo";$("start").value=10001;$("end").value=12400;$("motif").value="CG";$("flank").value=100;$("step").value=1;$("motif-both").checked=false;$("guide-both").checked=true;
    $("guides").value=`gR1  ${seq.slice(360,380)}\ngR2  ${seq.slice(1090,1110)}\ngR3  ${C.reverseComplement(seq.slice(1810,1830))}`;
    $("bed-text").value="chrDemo\t10000\t10600\tregion_A\nchrDemo\t10600\t11600\tregion_B\nchrDemo\t11600\t12400\tregion_C";
    pasteLength();return run(async()=>{await analyze();await markGuides();notice("サンプル解析完了 · CpGの分布と3本のgRNAを表示しています。解析方法をBED領域に切り替えて比較できます。");});
  }
  document.querySelectorAll("[data-source],[data-refmode],[data-mode]").forEach(b=>b.addEventListener("click",()=>{const key=Object.keys(b.dataset)[0];if(state[key]===b.dataset[key])return;state[key]=b.dataset[key];dirty(key!=="mode");sync();}));
  $("analyze").addEventListener("click",()=>run(analyze));$("mark-guides").addEventListener("click",()=>run(markGuides));
  $("bw-label").addEventListener("input",()=>{
    try{localStorage.setItem("motifStudio.bwLabel",$("bw-label").value);}catch{}
    if(state.result&&!state.stale)refreshFigureLabels(state.result);
    if(!state.stale)renderChart();
  });
  $("genome-version").addEventListener("change",()=>run(async()=>{
    if(state.bw?.genome&&state.bw.genome!==$("genome-version").value){state.bw=null;$("bw-info").textContent="ゲノム版が変わったためbigWigを読み込み直してください。";$("chip-atlas-status").textContent="選択中のゲノム版に対応するChIP/ATAC/DNase-seqのbigWigを取得します。";}
    state.genome=null;state.reference=null;state.location=null;state.locationCandidates=[];
    $("location-results").hidden=true;$("genome-file").value="";$("ref-file").value="";
    $("genome-info").textContent="";$("ref-info").textContent="ファイル未選択";
    $("location-status").textContent="参照ゲノムが変更されました。もう一度検索してください。";
    dirty(true);sync();
    if($("genome-version").value==="custom"){$("genome-status").textContent="検索用または領域解析用の参照FASTAを指定してください。";notice("自分のFASTAを使用します。");}
    else {await ensureGenome();notice("参照ゲノムを選択しました。配列の自動検索または領域解析を実行できます。");}
  }));
  $("locate").addEventListener("click",()=>run(locate));
  $("apply-location").addEventListener("click",()=>run(async()=>{const match=state.locationCandidates[Number($("location-match").value)];if(!match)throw new Error("検索をやり直してください。");await applyLocation(match,match.query);}));
  $("query-file").addEventListener("change",()=>run(async()=>{const file=$("query-file").files[0];if(!file)return;if(file.size>5*1024**2)throw new Error("検索配列ファイルは5 MB以下にしてください。");const text=await file.text();C.parseSequence(text);$("sequence").value=text;state.location=null;dirty(true);pasteLength();notice("検索する配列を読み込みました。参照ゲノムを指定して自動検索できます。");}));
  $("genome-file").addEventListener("change",()=>run(async()=>{const file=$("genome-file").files[0];if(!file)return;state.genome=null;state.location=null;state.locationCandidates=[];$("location-results").hidden=true;dirty(true);const data=await upload("/api/reference",file);state.genome={...data,name:file.name};$("genome-info").textContent=`${file.name} · ${Object.keys(data.chroms).length}配列`;notice("検索用参照を読み込みました。「配列を自動検索」で座標を探します。");}));
  $("demo").addEventListener("click",demo);$("empty-demo").addEventListener("click",demo);
  $("help-toggle").addEventListener("click",()=>{$("help").hidden=!$("help").hidden;$("help-toggle").setAttribute("aria-expanded",String(!$("help").hidden));});
  $("sequence").addEventListener("input",()=>{dirty(true);pasteLength();});
  ["chrom","start","end"].forEach(id=>$(id).addEventListener("input",()=>{dirty(true);if(id==="start")pasteLength();}));
  ["motif","motif-both","flank","step","bin-size","bed-text"].forEach(id=>$(id).addEventListener("input",()=>{dirty();if(id==="motif"||id==="motif-both"){state.scan=null;state.scanMotif=null;if(state.seq)renderSequence();}sync();}));
  ["bw-independent","bw-mode","bw-flank","bw-step","bw-bin-size","bw-bed-text"].forEach(id=>$(id).addEventListener("input",()=>{dirty();sync();}));
  $("bw-bed-file").addEventListener("change",()=>run(async()=>{const file=$("bw-bed-file").files[0];if(!file)return;if(file.size>20*1024**2)throw new Error("BEDは20 MB以下にしてください。");$("bw-bed-text").value=await file.text();dirty();C.parseBed($("bw-bed-text").value);notice("bigWig用BEDを読み込みました。「解析する」で更新してください。");}));
  for(const track of ["motif","bigwig"])$("table-"+track).addEventListener("click",()=>{state.tableTrack=track;state.page=0;renderTable();});
  ["guides","guide-both"].forEach(id=>$(id).addEventListener("input",()=>{state.guideResult=null;$("guide-results").hidden=true;$("metric-four").textContent="—";$("guide-status").textContent="gRNAが変更されています。「決定・マーク」で更新してください。";if(state.seq){renderSequence();renderChart();}buttons();}));
  $("cpg-preset").addEventListener("click",()=>{$("motif").value="CG";state.scan=null;state.scanMotif=null;dirty();if(state.seq)renderSequence();});
  $("bed-file").addEventListener("change",()=>run(async()=>{const file=$("bed-file").files[0];if(!file)return;if(file.size>20*1024**2)throw new Error("BEDは20 MB以下にしてください。");$("bed-text").value=await file.text();dirty();C.parseBed($("bed-text").value);notice(`${file.name} を読み込みました。「解析する」で集計してください。`);}));
  $("ref-file").addEventListener("change",()=>run(async()=>{const file=$("ref-file").files[0];if(!file)return;dirty(true);state.reference=null;const data=await upload("/api/reference",file);state.reference={...data,name:file.name};const [chrom,length]=Object.entries(data.chroms)[0];$("chrom").value=chrom;$("start").value=1;$("end").value=Math.min(length,10000);$("ref-info").textContent=`${file.name} · ${Object.keys(data.chroms).length}配列 · ${fmt(Object.values(data.chroms).reduce((a,b)=>a+b,0))} bp`;sync();notice("参照FASTAを読み込みました。表示領域を指定して解析してください。BED集計は全参照を対象にできます。");}));
  $("chip-atlas-fetch").addEventListener("click",()=>run(importChipAtlas));
  $("chip-atlas-id").addEventListener("keydown",event=>{if(event.key==="Enter"){event.preventDefault();run(importChipAtlas);}});
  $("bw-file").addEventListener("change",()=>run(async()=>{const file=$("bw-file").files[0];if(!file)return;state.bw=null;const data=await upload("/api/bigwig",file);useBigWig(data,file.name);$("chip-atlas-status").textContent="選択中のゲノム版に対応するChIP/ATAC/DNase-seqのbigWigを取得します。";notice("bigWigを読み込みました。解析方法と領域を選んで解析してください。");}));
  $("save-csv").addEventListener("click",()=>exportResults("csv"));$("save-bedgraph").addEventListener("click",()=>exportResults("bedgraph"));$("save-bed").addEventListener("click",()=>exportResults("bed"));
  $("save-bw-bedgraph").addEventListener("click",()=>exportResults("bw-bedgraph"));
  $("save-figure").addEventListener("click",()=>run(saveFigure));
  $("save-guides").addEventListener("click",()=>{if(state.guideResult&&state.seq)download(state.guideResult.hits.map(h=>[state.seq.chrom,h.start,h.end,safeName(h.name),0,h.strand].join("\t")).join("\n")+"\n","gRNA_matches.bed");});
  $("table-prev").addEventListener("click",()=>{state.page--;renderTable();});$("table-next").addEventListener("click",()=>{state.page++;renderTable();});
  $("seq-prev").addEventListener("click",()=>{state.seqPage--;renderSequence();});$("seq-next").addEventListener("click",()=>{state.seqPage++;renderSequence();});
  try{$("bw-label").value=localStorage.getItem("motifStudio.bwLabel")||"bigWig";}catch{}
  sync();
  if(location.protocol==="http:"||location.protocol==="https:")fetch("/api/health").then(r=>r.json()).then(data=>{state.server=true;$("session-info").hidden=false;restoreStartup(data.startup);genomeCatalog().then(entries=>{const entry=entries.find(e=>e.assembly===$("genome-version").value);if(entry)$("genome-status").textContent=genomeStatus(entry);}).catch(e=>notice(e.message,"warn"));$("connection").textContent=data.bigwig?"● ローカル接続":"● 配列解析のみ";$("connection").classList.add("connected");if(!data.bigwig)notice("bigWigを使用する場合は、サーバーの環境に requirements.txt の依存関係をインストールしてください。","warn");}).catch(()=>offline());else offline();
  function offline(){$("connection").textContent="HTMLモード";notice("配列の貼り付け解析を使用できます。全ゲノムFASTA・bigWigを扱う場合はローカルサーバーを起動してください。","warn");}
  function restoreStartup(data){
    if(!data||$("sequence").value.trim())return;
    if(data.genome){$("genome-version").value="custom";$("genome-status").textContent="起動時に指定した参照FASTAを使用します。";state.genome=data.genome;$("genome-info").textContent=data.genome.name;}
    if(data.bw){state.bw=data.bw;$("bw-info").textContent=`${data.bw.name} · ${Object.keys(data.bw.chroms).length}染色体`;state.source="overlay";}
    if(data.query)$("sequence").value=data.query;
    pasteLength();
    if(data.genome&&data.query)run(locate);
  }
})();
