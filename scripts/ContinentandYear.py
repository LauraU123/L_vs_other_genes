#!/usr/bin/python

import pandas as pd
import csv
from Bio import SeqIO
import matplotlib as plt

def GetContinent(country):
    asia = ['South Korea','Afghanistan', 'China','Saudi_Arabia', 'Kuwait', 'Qatar','Philippines', 'Jordan','Japan','Thailand','Vietnam','India','South_Korea',
            'Iran','Taiwan','Nepal','Lebanon','Hong_Kong','Myanmar','Mongolia','Pakistan','Malaysia','Cambodia','Israel','Laos']
    europe = ['Germany','Spain', 'France', 'Italy', 'Netherlands', 'Norway', 'Sweden','Czech Republic', 'Finland',
      'Denmark', 'United_Kingdom', 'Czech Republic', 'Switzerland', 'UK', 'England', 'Poland', 'Greece','Austria',
      'Bulgaria', 'Hungary', 'Luxembourg', 'Romania' , 'Slovakia', 'Estonia', 'Slovenia','Portugal',
      'Croatia', 'Lithuania', 'Latvia','Serbia', 'Estonia', 'ME', 'Iceland','Turkey','Russia','Belgium','Ireland']
    africa = ['Kenya',"Cote_d'Ivoire",'Central_African_Republic','Morocco','Cameroon', 'Nigeria', 'Tunisia', 'Africa','Madagascar', 'Gambia', 'Tanzania','Zambia','South_Africa','Ghana','Sudan','Mali','Egypt','Mozambique']
    north_america = ['USA', 'Canada','Mexico']
    south_america =['Panama','Brazil','Argentina','Guatemala','Cuba','Peru','Paraguay','Nicaragua','Colombia','Ecuador','Uruguay','Bolivia']
    oceania = ['Australia','New_Zealand','Fiji','NewZealand']
    if country in asia:
        return "Asia"
    elif country in europe:
        return "Europe"
    elif country in north_america:
        return "North America"
    elif country in south_america:
        return "South America"
    elif country in africa:
        return "Africa"
    elif country in oceania:
        return "Oceania"
    else:
        return "other"

with open ('longsequencesaligned.fasta') as f, open ('/home/laura/rsvforsnakemake/data/'+'metadata.tsv') as g:
    rd = csv.reader(g, delimiter="\t", quotechar='"')
    Continents =('Europe', 'Asia', 'Africa', 'North America', 'South America', 'Oceania')
    colours = ('skyblue', 'violet', 'forestgreen', 'indigo', 'purple','deepskyblue')
    file = SeqIO.parse(f, 'fasta')
    listofgenes=[]
    listofcontinents=[]
    listoflists=[]
    listofcountries=[]
    for record in file:  listofgenes.append(record.id)
    for row in rd:

        a = (row[-1].split('-'))
        b =row[0]

        if b in listofgenes:
            listofcountries.append(row[3])
            listoflists.append(row[-1])

    for country in listofcountries: listofcontinents.append(GetContinent(country))

    df = pd.DataFrame({'Date':listoflists, 'Continent':listofcontinents})
    grouped_df = df.groupby(df.Continent)
    df_dictionary=dict()

    for continent, colour in zip(Continents, colours):

        a = grouped_df.get_group(continent)
        del a['Continent']
        a["Date"] = a["Date"].astype("datetime64")
        fig1 = a.groupby(a["Date"].dt.year).count().plot(kind="bar", title='Sequences per Year '+str(continent)+'', xlabel ='year', color=colour)
        fig2 = a.groupby([a["Date"].dt.year, a["Date"].dt.month]).count().plot(kind="bar",figsize=(20, 5), color=colour, xlabel='year, month', title='Sequences per Month '+str(continent)+'')
        df_dictionary[str(continent)] = a
        fig1.figure.savefig("GraphYear" + str(continent) +".png", format="PNG")
        fig2.figure.savefig("GraphMonth" + str(continent) +".png", format="PNG")
