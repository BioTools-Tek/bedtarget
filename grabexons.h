#ifndef GRABEXONS_H
#define GRABEXONS_H

#include <QProcess>
#include <qstring.h>
#include <QList>
#include <QMap>
#include <iostream>

#include "splice.h"

// Gene1|Exon5_5'UTR

#define Exn "Exon"
#define Itr "Intron"
#define Ing "Intergenic"
#define AccSpl "SpliceA"
#define DonSpl "SpliceD"
#define UTR5 "5'UTR"
#define UTR3 "3'UTR"



using namespace std;

class ExonHolder{
//Chrom, genename, and direction are held by the contained GeneHolder
public:
    ExonHolder(uint cds1, uint cds2, uint exn, short ut, short fram){
        start = cds1; stop = cds2; exon = exn; utr=ut; frame=fram;
    }

    uint start, stop, exon;
    short utr;// -1 = not, 3 = 3'UTR, 5= 5'UTR

    //Exon Frames give an indication of how one exon joins into the other
    //e.g. 0 = reading frame starts at exon, 1 = reading frame carries 1 over from previous exon, 2 = carries 2 over from previous exon
    short frame;
};

class GeneHolder{

public:

    GeneHolder(QString name){gene_name = name;}
    QString gene_name;
    bool direction;  // '+' '-'
    QString chrom;

    int exonCount; //Needs to agree with size of start and stop

    //Defines coding regions
    uint cdsStart;
    uint cdsStop;

    uint txStart;
    uint txStop;
    QChar score_start;
    QChar score_endl;

    QList<ExonHolder *> exons;
    QList<uint>  exon_numbers;

    short reading_frames[500]; // Genes are very unlikely to have 500 exons, enough right?
    // http://biology.stackexchange.com/questions/29/what-are-the-limiting-factors-for-gene-length-and-number-of-exons

    uint determineExonNumber(uint exon){
        if (direction) return exon;
        else {
            int index_of = exon_numbers.indexOf(exon);

            //compute new index counting from the back of list
            index_of = (exon_numbers.length()-index_of) - 1;

            return exon_numbers.at(index_of);
        }
        cerr << "Could not determine direction!" << endl;
        exit(-1);
    }
};


class GrabExons{
public:
    QString database;
    QString table;
    bool local;
    QString chrom;
    int reg1, reg2;

    QMap<QString, QMap<QString, uint> > chromemap;
    QList<GeneHolder*> genes;
    //Map[Chrome][Gene][Exon] = #No. times this exon has been read already.

    //Chain
    GrabExons(QString database, QString region) {GrabExons(database, region, QString("refGene"));}
    GrabExons(QString database, QString region, QString table){ GrabExons(database, region, table, false);}
    GrabExons(QString database, QString region, QString table, bool local){
        this->database = database;
        this->table = table;
        this->local = local;

        QStringList reg = region.split(':');
        //if (reg.at(0).split("hr").at(1).toInt()==0){ cerr << "Invalid chromosome given" << endl; exit(-1);}
        this->chrom = reg.at(0).trimmed();

        this->reg1 = -1; this->reg2 = -1;

        if (reg.length()==2){
            QStringList posit = reg.at(1).split('-');
            if(posit.length()!=2){ cerr << "Please specify a min and max region" << endl; exit(-1);}
            this->reg1 = posit.at(0).toInt();
            this->reg2 = posit.at(1).toInt();
        }

        QList<QString> data = retrieveFromMYSQL();
        for (int i=0; i < data.size(); i++)
        {
            const QString &line = data.at(i);
            if(line.length()<5) continue;

            GeneHolder *toast = parseExons(line);
            if (toast->chrom!="0") genes.append(toast);
        }
    }

    QList<QString> retrieveFromMYSQL();
    uint isoformcheck(GeneHolder *gh);
    GeneHolder * parseExons(const QString &lines);
};


#endif // GRABEXONS_H
