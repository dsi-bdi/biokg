FROM python:3-slim

COPY . /biokg
WORKDIR /biokg
RUN pip install --no-cache-dir -r requirements.txt

CMD ["python", "run_all.py"]
