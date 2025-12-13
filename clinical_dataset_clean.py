import pandas as pd
import numpy as np
import re
df = pd.read_csv("CML Clinical DB(Variables Clínicas pacientes LM).csv", encoding = "latin1", sep = ';')

# Selecció i renom de columnes
pacients = df[
    ["ID",
     "RM profunda (BCR/ABL < 001%)",
     "RMM (BCRABL < 01%)",
     "BCR/ABL mes 3 (%)"]
].rename(
    columns={
        "RM profunda (BCR/ABL < 001%)": "RMP",
        "RMM (BCRABL < 01%)": "RMM",
        "BCR/ABL mes 3 (%)": "RM3"
    }
)

# Conversió segura RM3 → float
pacients["RM3"] = (
    pacients["RM3"]
    .astype(str)
    .str.strip()
    .str.replace(",", ".", regex=False)
    .replace(r"^NA.*$", np.nan, regex=True)
    .replace(["nan", "None", ""], np.nan)
)

pacients["RM3"] = pd.to_numeric(pacients["RM3"], errors="coerce")

# Subgrups clínics
pacients_rmp = df[pacients["RMP"] == 1]

pacients_rmm = df[
    (pacients["RMM"] == 1) &
    (pacients["RMP"] == 0)
]

pacients_rm3 = df[
    (pacients["RMM"] == 0) &
    (pacients["RMP"] == 0) &
    (pacients["RM3"] <= 10)
]

pacients_no_rm = df[
    (pacients["RMM"] == 0) &
    (pacients["RMP"] == 0) &
    (pacients["RM3"] > 10)
]

#Filtrar per Recaiguda ----------------------------------------------------


# Columnes de BCR/ABL
bcr_cols = [
    col for col in df.columns
    if col.startswith("BCR/ABL mes")
]

# Ordenar columnes pels mesos (3, 6, 12, 18, 24, 60)
def get_month(col):
    return int(re.search(r"mes (\d+)", col).group(1))

bcr_cols = sorted(bcr_cols, key=get_month)

# 1️ Normalitzar valors: strings → float, , → .
for col in bcr_cols:
    df[col] = (
    df[col]
    .astype(str)
    .str.strip()
    .str.replace(",", ".", regex=False)
    .replace(r"^NA.*$", np.nan, regex=True)
    .replace(["nan", "None", ""], np.nan)
    .astype(float)
)


# 2️ Detectar pujades
def hi_ha_pujada(row):
    valors = row[bcr_cols].values
    for i in range(len(valors) - 1):
        if (
            not np.isnan(valors[i])
            and not np.isnan(valors[i + 1])
            and valors[i + 1] > valors[i]
        ):
            return True
    return False

df["pujada_BCR_ABL"] = df.apply(hi_ha_pujada, axis=1)

# 3️ Opcional: detectar entre quins mesos
def mesos_pujada(row):
    pujades = []
    valors = row[bcr_cols].values
    for i in range(len(valors) - 1):
        if (
            not np.isnan(valors[i])
            and not np.isnan(valors[i + 1])
            and valors[i + 1] > valors[i]
        ):
            m1 = get_month(bcr_cols[i])
            m2 = get_month(bcr_cols[i + 1])
            pujades.append(f"{m1}→{m2}")
    return ", ".join(pujades)

df["mesos_pujada"] = df.apply(mesos_pujada, axis=1)

# 4️ Resultat final

recaiguda = df[df["pujada_BCR_ABL"] == True]
no_recaiguda = df[df["pujada_BCR_ABL"] == False]


# Crear còpia de la database "neta" --------------------------------------------------------------

import pandas as pd
import sqlite3

# Crear còpia del df original
df2 = df.copy()

# Inicialitzar columnes noves
df2["resposta_molecular"] = pd.NA
df2["recaiguda"] = pd.NA

# Assignar resposta molecular
df2.loc[df2["ID"].isin(pacients_rmp["ID"]), "resposta_molecular"] = "RMP"
df2.loc[df2["ID"].isin(pacients_rmm["ID"]), "resposta_molecular"] = "RMM"
df2.loc[df2["ID"].isin(pacients_rm3["ID"]), "resposta_molecular"] = "RM3"
df2.loc[df2["ID"].isin(pacients_no_rm["ID"]), "resposta_molecular"] = "No RM"

# Assignar recaiguda
df2.loc[df2["ID"].isin(recaiguda["ID"]), "recaiguda"] = "sí"
df2.loc[df2["ID"].isin(no_recaiguda["ID"]), "recaiguda"] = "no"



# Més neteja (formatejar dates, calcular edat i grups d'edat, ...) ---------------

from datetime import date

# Convertir dates
def convertir_data(x):
    if pd.isna(x):
        return pd.NaT

    x = str(x).strip().lower()

    try:
        dia, mes_txt, any_txt = x.split("-")
        dia = int(dia)
        any = int(any_txt)

        # Any 2 dígits → 4 dígits
        if any > 25:
            any = 1900 + any
        else:
            any = 2000 + any

        mesos = {
            "ene": 1, "feb": 2, "mar": 3, "abr": 4,
            "may": 5, "jun": 6, "jul": 7, "ago": 8,
            "sep": 9, "oct": 10, "nov": 11, "dic": 12
        }

        mes = mesos.get(mes_txt)

        if mes is None:
            return pd.NaT

        return pd.Timestamp(year=any, month=mes, day=dia)

    except Exception:
        return pd.NaT

df2["Data naixement"] = df2["Data naixement"].apply(convertir_data)
df2["Data naixement"] = pd.to_datetime(df2["Data naixement"], errors="coerce")
today = pd.Timestamp(date.today())

#Crear columnes per l'edat i el grup segons l'edat
df2["edat"] = today.year - df2["Data naixement"].dt.year


bins = [0, 40, 60, 80, np.inf]
labels = ["<40", "40–59", "60-79", "≥80"]

df2["edat_grup"] = pd.cut(
    df2["edat"],
    bins=bins,
    labels=labels,
    right=False
)


#Classificar segons el sexe
df2["SEXO"] = df["SEXO"].map({
    0: "Dona",
    1: "Home"
})

#Reformatejar altres dates incorrectes

df2["Fecha primera muestra"] = df2["Fecha primera muestra"].apply(convertir_data)
df2["Data Diagnòstic"] = df2["Data Diagnòstic"].apply(convertir_data)
df2["Data Alo-TPH"] = df2["Data Alo-TPH"].apply(convertir_data)


def convertir_data_ddmmaaany(x):
    try:
        # Passar a string si és int
        s = str(x)

        # Si el format és 'ddmmaaany' sense separadors, afegim '/'
        if s.isdigit() and len(s) == 8:
            s = s[:2] + '/' + s[2:4] + '/' + s[4:]

        # Convertir a datetime, format dia/mes/any
        return pd.to_datetime(s, format='%d/%m/%Y', errors='coerce')
    except:
        return pd.NaT

df2["INICIO ITC 1Âª LINEA"] = df2["INICIO ITC 1Âª LINEA"].apply(convertir_data_ddmmaaany)
df2["FIN ITC 1Âª LINEA"] = df2["FIN ITC 1Âª LINEA"].apply(convertir_data_ddmmaaany)
df2["INICIO ITC 2Âª LINEA"] = df2["INICIO ITC 2Âª LINEA"].apply(convertir_data_ddmmaaany)
df2["FIN ITC 2Âª LINEA"] = df2["FIN ITC 2Âª LINEA"].apply(convertir_data_ddmmaaany)
df2["INICIO ITC 3Âª LINEA"] = df2["INICIO ITC 3Âª LINEA"].apply(convertir_data_ddmmaaany)
df2["FIN ITC 3Âª LINEA"] = df2["FIN ITC 3Âª LINEA"].apply(convertir_data_ddmmaaany)
df2["INICIO ITC 4Âª LINEA"] = df2["INICIO ITC 4Âª LINEA"].apply(convertir_data_ddmmaaany)
df2["Fecha RMM"] = df2["Fecha RMM"].apply(convertir_data_ddmmaaany)
df2["Fecha RM profunda"] = df2["Fecha RM profunda"].apply(convertir_data_ddmmaaany)
df2["Fecha Stop ITC x RM profunda"] = df2["Fecha Stop ITC x RM profunda"].apply(convertir_data_ddmmaaany)
df2["Fecha Ultima Visita"] = df2["Fecha Ultima Visita"].apply(convertir_data_ddmmaaany)

# Canviar i uniformitzar el nom de les columnes
df3 = df2[
    ["Hospital", 
     "ID", 
     "FASTQ OK", 
     "Fecha primera muestra", 
     "Data Diagnòstic", 
     "Data naixement", 
     "SEXO", 
     "TIPO TRANSCRITO", 
     "Data Alo-TPH", 
     "ITC DE 1Âª LINEA", 
     "INICIO ITC 1Âª LINEA",
     "FIN ITC 1Âª LINEA", 
     "Causa Fin 1Âº ITC", 
     "ITC DE 2Âª LINEA", 
     "INICIO ITC 2Âª LINEA",
     "FIN ITC 2Âª LINEA",
     "Causa Fin 2Âª LINEA", 
     "ITC DE 3Âª LINEA", 
     "INICIO ITC 3Âª LINEA",
     "FIN ITC 3Âª LINEA", 
     "Causa Fin 3Âº ITC", 
     "ITC DE 4Âª LINEA", 
     "INICIO ITC 4Âª LINEA",
     "BCR/ABL al diagnÃ³stico (%)", 
     "BCR/ABL mes 3 (%)", 
     "BCR/ABL mes 6 (%)", 
     "BCR/ABL mes 12 (%)",
     "BCR/ABL mes 18 (%)", 
     "BCR/ABL mes 24 (%)", 
     "BCR/ABL mes 60 (%)", 
     "RMM (BCRABL < 01%)", 
     "Fecha RMM",
     "RM profunda (BCR/ABL < 001%)", 
     "Fecha RM profunda", 
     "Fecha Stop ITC x RM profunda",
     "Fecha Ultima Visita", 
     "Estat_UltimaVisita", 
     "Exitus", 
     "pujada_BCR_ABL", 
     "mesos_pujada",
     "resposta_molecular", 
     "recaiguda", 
     "edat", 
     "edat_grup"]    
].rename(
    columns={
        "Hospital" : "Hospital", 
        "ID" : "ID", 
        "FASTQ OK" : "FASTQ OK", 
        "Fecha primera muestra" : "Data primera mostra", 
        "Data Diagnòstic" : "Data diagnòstic", 
        "Data naixement" : "Data naixement ", 
        "SEXO" : "Sexe", 
        "TIPO TRANSCRITO" : "Tipus transcript", 
        "Data Alo-TPH" : "Data Alo-TPH", 
        "ITC DE 1Âª LINEA": "ITC 1a línia", 
        "INICIO ITC 1Âª LINEA" : "Inici ITC 1a línia",
        "FIN ITC 1Âª LINEA": "Fi ITC 1a línia", 
        "Causa Fin 1Âº ITC" : "Causa Fi ITC 1a línia", 
        "ITC DE 2Âª LINEA" : "ITC 2a línia", 
        "INICIO ITC 2Âª LINEA" : "Inici ITC 2a línia",
        "FIN ITC 2Âª LINEA" : "Fi ITC 2a línia",
        "Causa Fin 2Âª LINEA" : "Causa Fi ITC 2a línia", 
        "ITC DE 3Âª LINEA" : "ITC 3a línia", 
        "INICIO ITC 3Âª LINEA" : "Inici ITC 3a línia",
        "FIN ITC 3Âª LINEA": "Fi ITC 3a línia", 
        "Causa Fin 3Âº ITC": "Causa Fi ITC 3a línia", 
        "ITC DE 4Âª LINEA" : "ITC 4a línia", 
        "INICIO ITC 4Âª LINEA" : "Inici ITC 4a línia",
        "BCR/ABL al diagnÃ³stico (%)" : "BCR/ABL diagnòstic (%)", 
        "BCR/ABL mes 3 (%)" : "BCR/ABL mes 3 (%)", 
        "BCR/ABL mes 6 (%)" : "BCR/ABL mes 6 (%)", 
        "BCR/ABL mes 12 (%)" : "BCR/ABL mes 12 (%)",
        "BCR/ABL mes 18 (%)" : "BCR/ABL mes 18 (%)", 
        "BCR/ABL mes 24 (%)" : "BCR/ABL mes 24 (%)", 
        "BCR/ABL mes 60 (%)" : "BCR/ABL mes 60 (%)", 
        "RMM (BCRABL < 01%)" : "RMM (BCRABL < 01%)", 
        "Fecha RMM" : "Data RMM",
        "RM profunda (BCR/ABL < 001%)" : "RMP (BCR/ABL < 001%)", 
        "Fecha RM profunda" : "Data RMP", 
        "Fecha Stop ITC x RM profunda" : "Data Stop ITC x RMP",
        "Fecha Ultima Visita" : "Data Última Visita", 
        "Estat_UltimaVisita" : "Estat Última Visita", 
        "Exitus" : "Exitus", 
        "pujada_BCR_ABL" : "Pujada BCR/ABL", 
        "mesos_pujada" : "Mesos Pujada",
        "resposta_molecular": "Resposta Molecular", 
        "recaiguda" : "Recaiguda", 
        "edat": "Edat", 
        "edat_grup" : "Grup Edat"
    }
)

df3.to_csv(
    "CML_df_final.csv",
    index=False,
    encoding="utf-8"
)
