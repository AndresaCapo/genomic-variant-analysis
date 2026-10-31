FROM python:3.11-slim

WORKDIR /app

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

RUN pip install --no-cache-dir pysam==0.23.0

COPY analise_variantes.py .
COPY suporte_alelico.py .
COPY exportar_resultados.py .
COPY test_regressao.py .

ENTRYPOINT ["python3", "/app/exportar_resultados.py"]
