#ifndef PRINTEXONS_H
#define PRINTEXONS_H

#include "grabexons.h"
#include <string>

class PrintExons
{
public:
    PrintExons(QList<GeneHolder *> &genes, bool &exons, bool &introns, bool &utr, bool &intergenic, Splice *splice, bool &unique_names){
        uint last_genepos = 0; //Intergenic position holder

        for (int i=0; i < genes.size(); i++)
        {
            GeneHolder *gh = genes.at(i);

            uint last_exonpos = last_genepos;

            string chrom = gh->chrom.toStdString();
            string gene_name = gh->gene_name.toStdString();


            if (intergenic){
                cout << chrom << '\t'
                     << last_exonpos << '\t' << gh->txStart << '\t'
                     << gene_name << Ing << endl;
            }

            int exlength = gh->exons.length();


            for (int j=0; j < exlength; j++)
            {
                ExonHolder *ex = gh->exons.at(j);

                //#######  UTR -- SPLICE -- EXONS/UTR_Exon -- SPLICE #######//

                if(j!=0){   //Skip first
                    //first splice sites - Donor is 5', Acceptor is 3'
                    if(splice->splice_site!=0){
                        cout << chrom << '\t' << (ex->start - splice->splice_site)  << '\t' << ex->start << '\t'
                             << gene_name << Exn << ex->exon;

                        if(gh->direction && splice->acceptor_sites) cout << AccSpl << endl;
                        //reverse, Donor comes first instead
                        if(!gh->direction && splice->donor_sites) cout << DonSpl << endl;
                    }
                }

                ////// Exons || Exons_UTR
                if(exons){
                    cout << chrom << '\t'
                         << ex->start << '\t' << ex->stop << '\t'
                         << gene_name
                         << Exn << ex->exon << flush;
                }

                if(utr){
                    if(ex->utr!=-1){
                        if(ex->utr==5) cout << UTR5;
                        else cout << UTR3;
                    }
                }
                if(exons || utr) cout << endl;


                if(j!=(exlength-1)){   //Skip last
                    //second splice sites
                    if(splice->splice_site!=0){
                        cout << chrom << '\t' << ex->stop  << '\t' << (ex->stop + splice->splice_site) << '\t'
                             << gene_name << Exn << ex->exon << flush;

                        ex->stop += splice->splice_site;

                        if(gh->direction && splice->donor_sites) cout << DonSpl << endl;
                        //reverse, Acceptor comes last instead
                        if(!gh->direction && splice->acceptor_sites) cout << AccSpl << endl;
                    }
                }

                if(introns){
                    cout << chrom << '\t'
                         << last_exonpos << '\t' << (ex->start - splice->splice_site) << '\t'
                         << gene_name
                         << Itr << ex->exon << endl;
                }


                last_exonpos = ex->stop;
            }

            //Store last gene position
            last_genepos = gh->txStop;
        }
    }
};

#endif // PRINTEXONS_H
