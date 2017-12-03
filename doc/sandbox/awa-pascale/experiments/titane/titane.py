# -*- coding: utf-8 -*-
import pandas as pd
df = pd.read_excel("Ess4_bis_Titane_0.xlsx")
data = df[[ u"Déformation vraie", u"Contrainte vraie"]]
data.columns = ["strain", "stress"]
data.to_csv("../Ess4_bis_Titane_0.csv", index = False, sep = "\t")
