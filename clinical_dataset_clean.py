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

# Conversió segura RM3 → float RM3=Resposta Mol·lecular als 3 mesos
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

"""
pacients_rmp["KIR"] = pacients_rmp["ID"].map(
    df2.set_index("ID")["al·lels"]
)
pacients_rmm["KIR"] = pacients_rmm["ID"].map(
    df2.set_index("ID")["al·lels"]
)
pacients_rm3["KIR"] = pacients_rm3["ID"].map(
    df2.set_index("ID")["al·lels"]
)
pacients_no_rm["KIR"] = pacients_rm["ID"].map(
    df2.set_index("ID")["al·lels"]
)

"""

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


# 2️ Detectar pujades respecte els altres mesos
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

# 3️ Opcional: detectar entre quins mesos ha pujat aka recaiguda
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


# Crear CSV --------------------------------------------------------------

import pandas as pd
import sqlite3

# Crear còpia del df original per afegir noves columnes
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


df2.to_csv(
    "CML_df_final.csv",
    index=False,
    encoding="utf-8"
)


