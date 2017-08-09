#include "targeter.h"

#define SWAPPED "[SWAPPED]"


void Targeter::targetTranscriptStartSites(){

    MySQL::retrieveMYSQLResult("");
}


void Targeter::targetTFBS(){

}


void Targeter::targetPromoters(short margin_us, short margin_ds){
    QList<GeneHolder*> &genes = this->genelist;

    int gen_size = genes.size();
    for (int i=0; i< gen_size; i++)
    {
        GeneHolder *gh = genes.at(i);
        bool direct = gh->direction;

        // Promoters are found before transcription starts (in 5' or 3'), but UCSC does not list these
        // in a very consistent way:
        //    if -  txStart should be greater than txEnd
        //    if +  txStart should be less than txEnd
        // However what the tablebrowser seems to do for the most is to simply give the smaller upstream
        // value as the txStart, regardless of orientation.
        // Need to flag these.

        string chrom = gh->chrom.toStdString();
        string gene_name = gh->gene_name.toStdString();

        uint txSmall = gh->txStart;
        uint txBig  = gh->txStop;

        if (gh->txStop < gh->txStart){
            txBig   = gh->txStart;
            txSmall = gh->txStop;
        }


        //  Promoter definition
        //  <----------|---------->
        //    A promoter is just a region spanning the txSTART site.

        // 5' --> 3'
        if(direct){
            //Region before
            if (margin_us!=-1){
                uint txStart_before = txSmall - margin_us;
                cout << chrom << '\t'
                     << txStart_before << '\t' << txSmall << '\t'
                     << gene_name
                     << "|Promoter_" << margin_us << "upstream" << flush;
                printExtras(gh, -1);
            }

            //Region After
            if(margin_ds!=-1){
                uint txStop_after = txSmall + margin_ds;
                cout << chrom << '\t' << (txSmall + 1) << '\t' << txStop_after << '\t'
                     << gene_name
                     << "|Promoter_" << margin_ds << "dnstream" << flush;
                printExtras(gh, -1);

            }
        }

        // 3' --> 5'
        else {

            if (margin_us!=-1){  // upstream is to the right
                uint txRevStart = txBig + margin_us + 1;

                cout << chrom << '\t' << (txBig+1) << '\t' << txRevStart << '\t'
                     << gene_name
                     << "|Promoter_" << margin_us << "upstream" << flush;
                //<< (flag_q_txS?SWAPPED:"") << flush;
                printExtras(gh, -1);
            }

            if(margin_ds!=-1){  // downstream is to the left
                uint txRevStop = txBig - margin_ds;
                cout << chrom << '\t' << (txRevStop+1) << '\t' << txBig << '\t'
                     << gene_name
                     << "|Promoter_" << margin_ds << "dnstream" << flush;
                printExtras(gh, -1);
            }
        }
    }
}


void inline Targeter::printExtras(GeneHolder *&gh, short frame){
    gh->printDetails();

    // -99 is null
    if (this->frames){cout << '\t' << frame << flush;}
    cout << endl;
}



void Targeter::targetExons()
{
    QList<GeneHolder*> &genes = this->genelist;

    int gen_size = genes.size();
    for (int i=0; i< gen_size; i++)
    {
        GeneHolder *gh = genes.at(i);
        string chrom = gh->chrom.toStdString();
        string gene_name = gh->gene_name.toStdString();

        for(int j=0; j < gh->exons.length(); j++)
        {
            ExonHolder *ex = gh->exons.at(j);
            //Dont print UTR
            if (ex->utr==-1)
            {
                cout << chrom << '\t'
                     << ex->start << '\t' << ex->stop << '\t'
                     << gene_name << '|'
                     << Exn << gh->determineExonNumber(ex->exon) << flush; printExtras(gh, ex->frame);
            }
        }
    }
}


void Targeter::targetUTRExons()
{
    QList<GeneHolder*> &genes = this->genelist;

    int gen_size = genes.size();
    for (int i=0; i< gen_size; i++)
    {
        GeneHolder *gh = genes.at(i);
        string chrom = gh->chrom.toStdString();
        string gene_name = gh->gene_name.toStdString();

        for(int j=0; j < gh->exons.length(); j ++)
        {
            ExonHolder *ex = gh->exons.at(j);

            // Print UTR ONLY
            // No need to worry about splice sites

            if(ex->utr!=-1){

                cout << chrom << '\t'
                     << ex->start << '\t' << ex->stop << '\t'
                     << gene_name << '|'
                     << Exn << gh->determineExonNumber(ex->exon) << flush;

                if(ex->utr==5) cout << '_' << UTR5;
                else cout << '_' << UTR3;

                cout << flush; printExtras(gh, ex->frame);
            }
        }
    }
}


void Targeter::targetUTRRegions()
{
    QList<GeneHolder*> &genes = this->genelist;

    int gen_size = genes.size();
    for (int i=0; i< gen_size; i++)
    {
        GeneHolder *gh = genes.at(i);
        string chrom = gh->chrom.toStdString();
        string gene_name = gh->gene_name.toStdString();

        cout << chrom << '\t'
             << gh->txStart << '\t' << gh->cdsStart
             << '\t' << gene_name  << '|'
             << (gh->direction?UTR5:UTR3) << flush;printExtras(gh);

        cout << chrom << '\t'
             << gh->cdsStop << '\t' << gh->txStop
             << '\t' << gene_name  << '|'
             << (gh->direction?UTR5:UTR3) << flush; printExtras(gh);

    }
}



void Targeter::targetIntergenic(bool codingOnly)
{
    QList<GeneHolder*> &genes = this->genelist;
    uint last_genepos = 0; // Intergenic position holder
    QString last_gene = "";

    int gen_size = genes.size();
    for(int i=0; i < gen_size; i++)
    {
        GeneHolder *gh = genes.at(i);

        uint tmp_start = codingOnly?gh->cdsStart:gh->txStart;
        uint tmp_stop = codingOnly?gh->cdsStop:gh->txStop;

        if((i==0) || (last_gene.split("-ISO")[0] == gh->gene_name.split("-ISO")[0])){
            last_genepos = tmp_stop;
            last_gene = gh->gene_name;
            continue;
        }

        cout << gh->chrom.toUtf8().data() << '\t' << flush;

        //Isoforms can be fickle
        if(last_genepos > tmp_start){
            cout << tmp_start << '\t' << last_genepos
                 << '\t' << gh->gene_name.toUtf8().data() << '-' << last_gene.toUtf8().data() << flush;
        } else {
            cout << last_genepos << '\t' << tmp_start
                 << '\t' << last_gene.toUtf8().data() << '-' << gh->gene_name.toUtf8().data() << flush;
        }
        cout << '|' << Ing << flush; printExtras(gh);

        last_genepos = tmp_stop;
        last_gene = gh->gene_name;
    }
}



void Targeter::targetGenes(bool codingOnly)
{
    QList<GeneHolder*> &genes = this->genelist;

    int gen_size = genes.size();
    for(int i=0; i < gen_size; i++)
    {
        GeneHolder *gh = genes.at(i);

        uint tmp_start = codingOnly?gh->cdsStart:gh->txStart;
        uint tmp_stop = codingOnly?gh->cdsStop:gh->txStop;

        cout << gh->chrom.toUtf8().data() << '\t'
             << tmp_start << '\t' << tmp_stop  << '\t'
             << gh->gene_name.toUtf8().data()
             << flush; printExtras(gh);
    }
}


void Targeter::targetSpliceOnly(Splice *ss, bool no_utr)
{
    QList<GeneHolder*> &genes = this->genelist;

    int gen_size = genes.size();
    for (int i=0; i< gen_size; i++)
    {
        GeneHolder *gh = genes.at(i);
        string chrom = gh->chrom.toStdString();
        string gene_name = gh->gene_name.toStdString();


        for(int j=0; j < gh->exons.length(); j ++)
        {
            ExonHolder *ex = gh->exons.at(j);

            // If dont print utr, and exon is not utr, skip.. wat
//            if (no_utr && ex->utr!=-1) continue;

            if (ex->utr != -1) continue; // No splice sites on Exons that are only UTR


 //         if (ex->utr==-1)
 //         {

            //Print SPLICE ONLY
            //First exon, only splice on cdsStop side
            //Last exon, only splice on cdsStart side

            uint splice1_start = ex->start;
            uint splice1_end = ex->start;

            uint splice2_start = ex->stop;
            uint splice2_end = ex->stop;

            // Bedtarget recieves the exons in +'ve order
            // so the last exon is the same in any given direction

            // txStart cdsStart  first-exon
            //+       ______5'UTR OOOOOOOOOOO 3'E__A____5'S  --->
            //- <---  ______3'UTR OOOOOOOOOOO 5'E__D____3'S

            //               last exon     cdsStop    txStop
            //+         5'S__D___5'E OOOOOOOOOOO 3'UTR_____ --->
            //-  <---   3'E__A___3'S OOOOOOOOOOO 5'UTR_____

            //the placing of splice_sites is actually independent of direction.
            // it depends only on first and last index

            if (!(ex->start==gh->cdsStart)) splice1_start -= ss->splice_site;
            if (!(ex->stop==gh->cdsStop)) splice2_end += ss->splice_site;

            //Every now and then there are UTR_exons so we should ignore any
            //             txStart             cdsStart
            // +      5'UTR_______DonSP ooooooo|OOOOOOOOO AccSP_____
            //splice sites at cdsStart or cdsStop, but not any inbetween ones

            // The order of printing splice sites should still be splice1 followed
            // by splice2. Just the naming of donor and acceptor change.

            //if direction is +ve: donor_sites = splice1, accept_sites = splice2
            //if direction is -ve: accept_sites= splice1,  donor_sites = splice2

            QString utr_mention="";
            if(ex->utr==5) utr_mention = QString('_').append(UTR5);
            else if (ex->utr==3) utr_mention = QString('_').append(UTR3);




            if(gh->direction){
//                cerr << gene_name << "|" << Exn << gh->determineExonNumber(ex->exon) << "--" << ex->exon << "/" << j << endl;

                if(ss->donor_sites){
                    if (splice1_start != splice1_end){ // Only process non-UTR boundary splice sites

                        cout << chrom << '\t'
                             << splice1_start << '\t' << splice1_end << '\t'
                             << gene_name << '|'
                             << Exn << gh->determineExonNumber(ex->exon) << utr_mention.toUtf8().data() << '_' << DonSpl << flush; printExtras(gh, ex->frame);
                    }
                }
                if(ss->acceptor_sites){
                    if (splice2_start != splice2_end){ // Only process non-UTR boundary splice sites

                        cout << chrom << '\t'
                             << splice2_start << '\t' << splice2_end << '\t'
                             << gene_name << '|'
                             << Exn << gh->determineExonNumber(ex->exon) << utr_mention.toUtf8().data() << '_' << AccSpl << flush; printExtras(gh, ex->frame);

                    }
                }
            }
            else{
//                cerr << gene_name << "|" << Exn << gh->determineExonNumber(ex->exon) <<  utr_mention.toUtf8().data() << '_' << "--" << ex->exon << "/" << j << endl;

                if(ss->acceptor_sites){
                    if (splice2_start != splice2_end){ // Only process non-UTR boundary splice sites

                        cout << chrom << '\t'
                             << splice2_start << '\t' << splice2_end << '\t'
                             << gene_name << '|'
                             << Exn << gh->determineExonNumber(ex->exon) << utr_mention.toUtf8().data() << '_' << AccSpl << flush; printExtras(gh, ex->frame);
                    }
                }
                if(ss->donor_sites){
                    if (splice1_start != splice1_end){ // Only process non-UTR boundary splice sites

                        cout << chrom << '\t'
                             << splice1_start << '\t' << splice1_end << '\t'
                             << gene_name << '|'
                             << Exn << gh->determineExonNumber(ex->exon) << utr_mention.toUtf8().data() << '_' << DonSpl << flush; printExtras(gh, ex->frame);
                    }
                }
            }
        }
    }
}


void Targeter::targetIntrons(Splice *ss)
{
    QList<GeneHolder*> &genes = this->genelist;

    int gen_size = genes.size();
    for (int i=0; i< gen_size; i++)
    {
        GeneHolder *gh = genes.at(i);
        string chrom = gh->chrom.toStdString();
        string gene_name = gh->gene_name.toStdString();

        for(int j=0; j < gh->exons.length(); j ++)
        {
            ExonHolder *ex = gh->exons.at(j);

            // out of bounds left
            if (ex->stop <= gh->txStart){
                continue;
            }

            //out of bounds right
            if (ex->start >= gh->txStop){
                break;
            }

            // always skip first exon, only want end position
            if (j == 0) {
                continue;
            }

            uint left = gh->exons.at(j-1)->stop; // last exon pos
            uint right = ex->start;

            if (left == right){
                continue;
            }

            // If using splice sites in same map, then adjust margins
            if (ss->donor_sites){    left  += ss->splice_site;}
            if (ss->acceptor_sites){ right -= ss->splice_site;}


            cout << chrom << '\t'
                 << left << '\t' << right << '\t'
                 << gene_name << '|'
                 << Itr << (gh->determineExonNumber(ex->exon) - 1) << flush; printExtras(gh, ex->frame);

        }
    }
}






