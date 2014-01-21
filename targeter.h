#ifndef TARGETER_H
#define TARGETER_H

#include "grabexons.h"
#include "Args.h"
#include <QList>


class Targeter
{
private:
    typedef struct {
        uint pos1;
        uint pos2;
    } pair;
public:
    QList<GeneHolder*> genelist;
    bool scores, frames, direction;

    QList<QChar> checks;

    Targeter(QList<GeneHolder*> &genelist, Args &args)
    {
        checks.append('n');
        checks.append('u');
        checks.append('i');
        checks.append('c');

        this->genelist = genelist;
        this->scores = args.scores;
        this->frames = args.exonframes;
        this->direction = args.direction;

        if(args.promoter_margin!=-1) targetPromoters(args.promoter_margin); // targets genes only

        if(args.intergenic) targetIntergenic(args.codingonly);

        if(args.introns) targetIntrons(args.ss);


        if(args.genes){

            if(args.utrlevel == 0) targetGenes(args.codingonly);
            else if(args.utrlevel == 1) {
                targetGenes(args.codingonly);
                targetUTRRegions();
            }
            else if(args.utrlevel == 2) targetUTRRegions();
        }

        if(args.exons){
            if(args.ss->splice_site!=0) {
                if (args.utrlevel == 1) targetSpliceOnly(args.ss,false); //no_utr=false
                else targetSpliceOnly(args.ss); //no_utr = true

            }

            if(!args.ss->splice_only){
                if(args.utrlevel == 0) targetExons();
                else if(args.utrlevel == 1){
                    targetExons();
                    targetUTRExons();
                }
            }
            else if(args.utrlevel == 2) targetUTRExons();

        }
    }

    //Target Jobs
    //Key:
    //    ~ -- coding in progress
    //    ' -- coding done, testing needed
    //    " -- testing
    //    ¬ -- bugs, needs fixing
    //    ^ -- working, complete

    //amount to subtract from middle exons
    //amount to subtract from middle exons
    void targetIntergenic(bool codingOnly);   //'
    void targetGenes(bool codingOnly);        //'
    void targetExons();                       //'
    void targetIntrons(Splice *ss);           //'
    void targetUTRExons();                    //'
    void targetUTRRegions();                  //'
    void targetSpliceOnly(Splice *ss, bool no_utr=true);        //'
    void targetSpliceUTR(Splice *ss);         //'
    void targetPromoters(short margin_ds);

    //Other jobs
    void sortedList();
    void printList();
    void printExtras(QChar &score1, QChar &score2, QChar myway=QChar('A'), short frame=-99);

};

#endif // TARGETER_H
