
# Predict to Protect! - Efecte de les variants al·lèliques de KIR en la recaiguda de pacients de Leucèmia Mieloide Crònica (CML)

Mireia Fàbregas, Pol Baulenas, Pepe Merino, Maria Pascual

Els objectius del repte són 
1. Dissenyar un software a través de la mostra (30 seqüències) per tal d’entendre com passar d’una seqüència a les variants al·lèliques.
2. Utilitzar el software per tal de deduir la variant al·lèlica de cada pacient
3. Deduir a través de la base de dades si hi ha una relació entre variants al·lèliques KIR:
   - La resposta molecular del tractament
   - La recaiguda del pacient
   - i la possibilitat de mantenir la remissió després de la retirada del tractament.
## 1. Software
## 2. Determinar variants al·lèliques
Hem utilitzat el KIRCLE per inferir els genotips dels gens KIR a partir de dades de seqüenciació d'alt rendiment (fitxers BAM). Per tant, hem hagut de fer un pas previ per preparar les dades i les bases de referència.
Malauradament, no hem pogut provar el software amb tota la cohort de pacients i controls per falta de temps. Són un gran nombre de dades i un software lent. 

### 2.1 Preparació de Dades i Bases de Referència
Hem creat dos shell scripts:
alinear.sh: Processa les dades de seqüenciació, convertint els fitxers fastq en un format d'alineament binari (BAM), necessari per a l'anàlisi posterior.
actualitzar_dbs.sh: És fonamental per a la precisió, ja que genera una base de dades al·lèlica de referència actualitzada i curada utilitzant el repositori oficial d'IPD-KIR (EMBL-EBI).

### 2.2 Ponderació Quadràtica per Resoldre l'Ambigüitat
L'algorisme EM (Expectation-Maximization) original del KIRCLE (2019) utilitzava un compte binari de lectures: cada alineació era assignada amb un pes d'1, independentment de la seva qualitat. Aquesta aproximació resultava problemàtica en gens altament polimòrfics i homòlegs com els KIR, ja que generava ambigüitat entre al·lels gairebé idèntics (que només difereixen en un SNP).

Per resoldre aquesta limitació, hem modificat la lògica central a EM_algorithm.py per substituir el compte binari per una ponderació quadràtica basada en el Bitscore de BLAST (\\(x=Bitscore^2\\)). Bitscore és una mètrica logarítmica que quantifica la qualitat d'una alineació, penalitzant mismatches i gaps. En elevar el Bitscore al quadrat, s'accentua de manera desproporcionada la diferència entre alineaments perfectes (puntuació alta) i alineaments amb un sol error (puntuació lleugerament inferior).

Aquest canvi converteix la qualitat de l'alineament en una mètrica de probabilitat molt més sensible. Penalitzem severament els al·lels que no coincideixen perfectament amb els reads, permetent a l'algorisme EM convergír més ràpidament i amb major seguretat cap al genotip correcte, distingint de forma efectiva entre variants al·lèliques gairebé idèntiques.

## 3. Anàlisi de la base de dades
### 3.1 Netejar dades
Per tal de deduir la relació entre les variants al·lèliques i el desenvolupament, hem de filtrar els pacients. Ho fem primer per resposta molecular (4 grups) i després per si han o no recaigut (2 grups). Hem utilitzat pandas, numpy i regular expressions, al fitxer clinical_dataset_clean.py.
### 3.2 Tria de pacients
Degut a la natura intensiva i requeridora del processament de dades bioinformàtiques, hem optat per seleccionar una subcohort representativa per a l'anàlisi, estratificada en quatre grups principals basats en la presència o absència de Resposta Molecular (RM). Dins de cada grup, hem seleccionat 8 pacients (un total de 32), mantenint una distribució equilibrada per gènere (4 dones i 4 homes), assegurant una variabilitat d'edats, i incloent una representació de l'evolució clínica (2 pacients amb recaiguda i 2 pacients sense recaiguda) per garantir la robustesa comparativa de la mostra.

### 3.3 Anàlisi del dataset i visualització de dades
Hem fet uns quants tests de Chi-squared per intentar determinar si hi ha alguna correlació entre algunes de les dades dels pacients, al fitxer correlation_tests.R. Per tal de visualitzar les dades, hem usat RStudio i ggplot, al fitxer data_visualization.R.

## 4. Model Predictiu
Hem desenvolupat un programa en C que processa dades al·lèliques KIR per identificar patrons de recaiguda en pacients amb CML. El software llegeix un arxiu de al·lels, extreu i agrupa les variants pel seu tipus principal, calculant la freqüència d'ocurrència i la probabilitat de recaiguda associada a cadascuna. Els al·lels es classifiquen ordenadament segons la seva probabilitat de recaiguda, generant dos arxius CSV: un amb tots els subtipos individuals i un altre amb els tipus principals consolidats. Això permet identificar quines variants al·lèliques de KIR presenten major risc de recaiguda i facilita la predicció del pronòstic clínic en futurs pacients amb CML basant-se en el seu genotip KIR.
