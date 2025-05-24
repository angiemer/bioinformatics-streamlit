# Χρήση επίσημης εικόνας Python
FROM python:3.10-slim

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

# Ορισμός working directory μέσα στο container
WORKDIR /app

# Εγκατάσταση system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libglib2.0-0 \
    libsm6 \
    libxext6 \
    libxrender-dev \
    && rm -rf /var/lib/apt/lists/*

# Αντιγραφή του αρχείου requirements.txt
COPY requirements.txt .

# Εγκατάσταση των Python dependencies
RUN pip install --upgrade pip
RUN pip install --no-cache-dir -r requirements.txt

# Αντιγραφή των αρχείων του project στον container
COPY . .

# Άνοιγμα της θύρας 8501 για την εφαρμογή Streamlit
EXPOSE 8501

# Εκκίνηση της εφαρμογής Streamlit
CMD ["streamlit", "run", "streamlit_app.py"]
