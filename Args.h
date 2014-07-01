#ifndef ARGS_H
#define ARGS_H

#define VERSION "1.9b"

#include <QString>
#include <iostream>
#include "splice.h"

#define assert(x){ cerr << x << endl; cerr << endl; exit(-1);}

using namespace std;

class Args
{
public:
    //DB Args
    QString table, database, region;
    bool local;

    //Gene Args
    bool genes, exons, codingonly;
    bool introns, intergenic, scores, exonframes, direction;
    short utrlevel, promoter_margin_downstream, promoter_margin_upstream;

    Splice *ss;

    Args(int argc, char * argv[])
    {
        //Defaults
        database = "hg19";
        table="refGene";
        local=false, codingonly=true;

        genes=false, exons=false;
        introns=false, intergenic = false;
        scores=false; exonframes=false;
        direction=false;

        promoter_margin_downstream = promoter_margin_upstream = -1;

        int splice_amount = 0;
        bool splice_don = true, splice_acc = true;
        bool utrtoo = false;
        bool utronly = false;
        bool spliceonly = false;

        if (argc > 1 ){
            if (QString(argv[1]).contains("chr")) { //region
                region = argv[1];
            } else {
                cerr << "Please specify a chromosome or region." << endl;
                usage();
            }

            //Opts
            bool errs=false;
            QString errorsT = "Cannot parse: ";
            QList<QString> okay_vars;

            for (int i=2; i < argc; i++){
                QString tmp = argv[i];

                //Map Options
                if(tmp=="--local") {local=true; okay_vars.append(tmp);}

                else if(tmp=="--table"){
                    if (argc <= i+1){
                        cerr << "Give a table name. Default is refGene." << endl;
                        usage();
                    }
                    table = argv[i+1];
                    okay_vars.append(argv[i+1]);
                    okay_vars.append(tmp);
                }

                else if(tmp=="--database"){
                    if (argc <= i+1){
                        cerr << "Give a database name. Default is hg19." << endl;
                        usage();
                    }
                    database = argv[i+1];
                    okay_vars.append(argv[i+1]);
                    okay_vars.append(tmp);
                }

                //Gene Options
                else if(tmp=="--genes") {genes=true; okay_vars.append(tmp);}
                else if(tmp=="--exons") {exons=true; okay_vars.append(tmp);}
                else if(tmp=="--include-noncoding") {codingonly=false; okay_vars.append(tmp);}
                else if(tmp=="--introns") {introns=true; okay_vars.append(tmp);}
                else if(tmp=="--intergenic") {intergenic=true; okay_vars.append(tmp);}
                else if(tmp=="--scores") {scores=true; okay_vars.append(tmp);}
                else if(tmp=="--frames") {exonframes=true; okay_vars.append(tmp);}
                else if(tmp=="--direction") {direction=true; okay_vars.append(tmp);}
                else if(tmp=="--utr") {utrtoo=true; okay_vars.append(tmp);}
                else if(tmp=="--utr-only") {utronly=true; okay_vars.append(tmp);}

                //Promoter margin opts Opts
                else if(tmp=="--promoters") {
                    if (argc <= i+2 ){
                        cerr << "Give downsream and upstream promoter values. Default is 5000 and 0." << endl;
                        exit(-1);
                    }
                    promoter_margin_downstream = QString(argv[i+1]).toShort();
                    promoter_margin_upstream = QString(argv[i+2]).toShort();
                    okay_vars.append(argv[i+1]);okay_vars.append(argv[i+2]);
                    okay_vars.append(tmp);
                }

                //Splice specific options
                else if(tmp=="--splice"){
                    if (argc <= i+1 || QString(argv[i+1]).toInt()==0 ){
                        cerr << "Give a splice value. Default is 0." << endl;
                        exit(-1);
                    }
                    splice_amount = QString(argv[i+1]).toInt();
                    okay_vars.append(argv[i+1]);
                    okay_vars.append(tmp);
                }
                else if(tmp=="--splice-only") {spliceonly=true; okay_vars.append(tmp);}
                else if(tmp=="--splice-no-acceptors") {splice_acc=false; okay_vars.append(tmp);}
                else if(tmp=="--splice-no-donors") {splice_don=false; okay_vars.append(tmp);}

                //Errors
                else {
                    bool error_yes=true;
                    for (int i=0; i< okay_vars.length(); i++){
                        if (tmp==okay_vars.at(i)){
                            error_yes=false;
                            break;
                        }
                    }
                    if (error_yes) {
                        errorsT.append(tmp+" ");
                        errs = true;
                    }
                }
            }
            if (errs) { cerr << errorsT.toUtf8().data() << endl; exit(-1);}
        }
        else usage();

        //Check opts
        if (genes && exons) assert("Either exons or genes, one or the other");

        if(!genes && !exons){ // Nothing specified
            if(!introns && !intergenic && promoter_margin_downstream==-1){
                usage();
//              assert("nothing specified");
            }
        }

        //Splice
        if(splice_amount!=0 && !exons) assert("--exons must be specified if --splice is used");
        if ((!splice_acc || !splice_don) && splice_amount==0) assert("--splice must be specified too");
        //UTR
        if (genes){
            if (utronly) utrlevel = 2;
            else if (utrtoo) utrlevel = 1;
            else utrlevel = 0;
        }
        if (exons){
            if (utronly) utrlevel = 2;
            else if (utrtoo) utrlevel = 1;
            else utrlevel = 0;
        }
        ss = new Splice(splice_amount,splice_don, splice_acc, spliceonly);
    }


    void usage(){
        cerr << "version " << VERSION << endl;
        cerr << "\nusage:" << "bedtarget2 " << "<chrN[:reg1-reg2]> [OPTIONS [MAP GENE]]" << endl;
        cerr << "\nMAP OPTIONS:" << endl;
        cerr << "   --database [hg18|hg19|etc]" << endl;
        cerr << "   --table [refGene]" << endl;
        cerr << "   --local " << endl;
        cerr << "\nGENE OPTIONS:" << endl;
        cerr << "   --utr" << endl;
        cerr << "   --utr-only" << endl;
        cerr << "   --promoters <downstream_margin> <upstream_margin>" << endl;
        cerr << "   --include-noncoding" << endl;
        cerr << "   --introns" << endl;
        cerr << "   --intergenic" << endl;
        cerr << "   --splice <margin> " << endl;
        cerr << "   --splice-only " << endl;
        cerr << "       --splice-no-acceptors" << endl;
        cerr << "       --splice-no-donors" << endl;
        cerr << "   --genes || --exons" << endl;
        cerr << "   --scores" << endl;
        cerr << "   --frames" << endl;
        cerr << "   --direction" << endl;
        cerr << endl;
        exit(-1);
    }

};



#endif // ARGS_H
