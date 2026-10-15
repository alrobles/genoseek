FROM python:3.12-slim

WORKDIR /app

COPY backend/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY backend/ .
COPY frontend/public/ /app/frontend/public/

EXPOSE 8651

CMD ["uvicorn", "api:app", "--host", "0.0.0.0", "--port", "8651"]
