# Bioinformatics Streamlit App

Διαδραστική εφαρμογή ανάλυσης δεδομένων Μοριακής Βιολογίας με Python και Streamlit.

## Περιγραφή

Η εφαρμογή παρέχει:

- Προεπεξεργασία δεδομένων scRNA-seq
- Ανάλυση PCA και clustering
- UMAP απεικονίσεις
- Εκπαίδευση μοντέλων Machine Learning
- Δυνατότητα παραμετροποίησης από τον χρήστη
- Tab με πληροφορίες για την ομάδα

## Εκκίνηση μέσω Docker

```bash
docker build -t my-streamlit-app .
docker run --rm -p 8501:8501 my-streamlit-app
