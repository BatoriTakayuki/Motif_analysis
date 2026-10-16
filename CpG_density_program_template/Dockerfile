FROM python:3.12-slim-bookworm
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY server.py index.html app.js core.js style.css ./
EXPOSE 8765
CMD ["python", "server.py", "--host", "0.0.0.0"]
