
# Predict to Protect! - Efecte de les variants al·lèliques de KIR en la recaiguda de pacients de Leucèmia Mieloide Crònica (CML)
Objectius
1. Dissenyar un software a través de la mostra (30 seqüències) per tal d’entendre com passar d’una seqüència a les variants al·lèliques.
2. Utilitzar el software per tal de deduir la variant al·lèlica de cada pacient
3. Deduir a través de la base de dades si hi ha una relació entre variants al·lèliques KIR:
   - La resposta molecular del tractament
   - La recaiguda del pacient
   - i la possibilitat de mantenir la remissió després de la retirada del tractament.
## 1. Software
## 2. Determinar variants al·lèliques
## 3. Anàlisi de la base de dades
### 3.1 Netejar dades
Per tal de deduir la relació entre les variants al·lèliques i el desenvolupament, hem de filtrar els pacients. Ho fem primer per resposta molecular (4 grups) i després per si han o no recaigut (2 grups). Hem utilitzat pandas, numpy i regular expressions, al fitxer clinical_dataset_clean.py.

### 3.2 Anàlisi del dataset i visualització de dades
Per tal de visualitzar les dades, hem usat RStudio i ggplot, al fitxer data_visualization.R.
