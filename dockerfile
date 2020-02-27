FROM python:3-alpine3.10

COPY . /biokg
WORKDIR /biokg
RUN pip install --no-cache-dir -r requirements.txt

CMD ["python", "main.py"]
