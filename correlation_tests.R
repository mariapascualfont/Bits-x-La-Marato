df <- read.csv("~/Documents/Clinical cases/CML_df_final.csv")

chisq.test(df$Resposta.Molecular, df$Recaiguda) #Chi-squared : 2.3983, p-value: 0.6629
# Això ens indica que no hi ha una relació directa entre la resposta molecular i la possibilitat de recaure,
# no perquè tenir una resposta molecular profunda sigui poc concloent (que ho és), sinó perquè tot i que recaiguis
# pots arribar a una resposta molecular profunda.

chisq.test(df$Resposta.Molecular, df$Tipus.transcript) #Chi-squared : 11.307, p-value: 0.5028
# La resposta molecular i el tipus de transcript no estan directament relacionats, hi ha algun altre factor
# que duu als diferents outcomes (possiblement les variants al·lèliques)

chisq.test(df$Resposta.Molecular, df$Grup.Edat) #Chi-squared : 18.106, p-value : 0.1125
chisq.test(df$Resposta.Molecular, df$Sexe) #Chi-squared : 1.02, p-value : 0.9062
chisq.test(df$Recaiguda, df$Grup.Edat) #Chi-squared : 0.58, p-value : 0.8996
chisq.test(df$Recaiguda, df$Sexe) #Chi-squared : 0, p-value : 1
chisq.test(df$Tipus.transcript, df$Grup.Edat) #Chi-squared : 10.191, p-value : 0.3352
chisq.test(df$Tipus.transcript, df$Sexe) #Chi-squared : 4.0874, p-value : 0.2522

# En conclusió, no hi ha cap relació aparent que suggereixi un biaix en la població o una covariant
# que jugui algun paper rellevant en la resposta molecular, la recaiguda dels pacients, o el 
# desenvolupament d'algun tipus de transcript més maligne. 

# L'únic factor que sembla una mica relacionat amb un altre és la Resposta Molecular segons el grup d'edat.
# Això és possible a que es degui a que 1. la població que estudiem és una població envellida, ja que els càncers
# apareixen més sovint en persones d'edat avançada, i 2. en que és probable que es dimensioni el tipus de 
# tractament a fer servir (si és més o menys agressiu), segons l'edat del pacient. 
