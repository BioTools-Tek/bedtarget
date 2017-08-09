#ifndef GRABEXONS_H
#define GRABEXONS_H

//#define DEBUG


#include <QProcess>
#include <qstring.h>
#include <QList>
#include <QMap>
#include <iostream>

#include "splice.h"
#include "mysql.h"

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

    static bool print_direction;  // '+' '-'
    static bool print_scores;
    static bool print_uniquenames;
    static bool print_frames;

    static void printHeader(){
        cout << "#CHROM\tSTART\tSTOP\tGENEINFO";
        if (print_scores) cout << "\tSCORE";
        if (print_direction) cout << "\tDIRECT";
        if (print_uniquenames) cout << "\tRGNAME";
        if (print_frames) cout << "\tFRAMES";
        cout << endl;
    }

    QString chrom;

    bool direction;

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

    QList<QChar> checks;
    QString gene_name;
    QString unique_name;

    short reading_frames[500]; // Genes are very unlikely to have 500 exons, enough right?
    // http://biology.stackexchange.com/questions/29/what-are-the-limiting-factors-for-gene-length-and-number-of-exons

    GeneHolder(QString name){
        checks.append('n');
        checks.append('u');
        checks.append('i');
        checks.append('c');
        gene_name = name;
    }

    uint determineExonNumber(uint exon)
    {
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


    void debuginfo(){
        cout << " [strand=" << (direction?'+':'-') << ",txStart=" << txStart
                                               << ",txStop=" << txStop << "]" << flush;
    }



    void printDetails(){
        if (print_scores){
            //Append Quality
            QList<QChar> &keys = checks;
            cout << '\t' << flush;

            if (score_start.isNull()){
                cout << "noscore" << flush;
            }
            else{
                if(keys.contains(score_start))  cout << score_start.toLatin1() << flush;
                else cout << 'w' << flush;

                if(keys.contains(score_endl))  cout << score_endl.toUpper().toLatin1() << flush;
                else cout << 'W' << flush;
            }
        }
        if (print_direction) cout << '\t' << (direction?'+':'-') << flush;
        if (print_uniquenames) cout << '\t' << unique_name.toUtf8().data() << flush;
    }
};


class GrabExons{
public:
    QString database;
    QString table;
    QString chrom;
    QMap<QString, int> position_map;

    bool local;
    int reg1, reg2;

    QMap<QString, QMap<QString, uint> > chromemap;
    QList<GeneHolder*> genes;
    //Map[Chrome][Gene][Exon] = #No. times this exon has been read already.

    //Chain
    GrabExons(QString database, QString region) {GrabExons(database, region, QString("refGene"));}
    GrabExons(QString database, QString region, QString table){ GrabExons(database, region, table, false);}
    GrabExons(QString database, QString region, QString table, bool local)
    {
        this->database = database;
        this->table = table;
        this->local = local;

        QStringList reg = region.split(':');
        //if (reg.at(0).split("hr").at(1).toInt()==0){ cerr << "Invalid chromosome given" << endl; exit(-1);}
        this->chrom = reg.at(0).trimmed();

        this->reg1 = -1; this->reg2 = -1;

        if (reg.length()==2){
            QStringList posit = reg.at(1).split('-');
//            if(posit.length()!=2){ cerr << "Please specify a min and max region" << endl; exit(-1);}

            bool rgb;
            this->reg1 = posit.at(0).toInt(&rgb);
            if (!rgb) cerr << "Please specify a start region for " << region.toUtf8().data() << endl; exit(-1);

            if (posit.size()==2){
                this->reg2 = posit.at(1).toInt(&rgb);
                if (!rgb) cerr << "Please specify a proper region for " << region.toUtf8().data() << endl; exit(-1);
            }
            else{ //assume end is max
                this->reg2 = 1 << 31;
            }
        }



        QList<QString> data = retrieveFromMYSQL();


        makePositionMap(data.at(0)); // Use header to get indexes

        for (int i=1; i < data.size(); i++)
        {
            const QString &line = data.at(i);
            if(line.trimmed().length()<5) continue;

            GeneHolder *toast = parseExons(line);
            if (toast->chrom!="0") genes.append(toast);

            // IDIOT! Don't delete the pointer, that's what you're storing! --> Destructor
            //delete toast;
        }
    }

    ~GrabExons(){
//        cerr << "Deleting pointers.. " << endl;
        for (GeneHolder *gh : genes) delete gh;
        genes.clear();
        chromemap.clear();
    }

    QList<QString> retrieveFromMYSQL();
    uint isoformcheck(GeneHolder *gh);
    GeneHolder * parseExons(const QString &lines);
    void makePositionMap(const QString &line);
};


#endif // GRABEXONS_H
