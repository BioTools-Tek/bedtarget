#ifndef SPLICE_H
#define SPLICE_H

class Splice {
public:
    int splice_site;
    bool donor_sites, acceptor_sites;
    bool splice_only;

    Splice(int splice_site){
        Splice(splice_site, true, true,false);
    }

    Splice(int splice_site, bool donor_sites, bool acceptor_sites, bool only){
        this->splice_site = splice_site;
        this->donor_sites = donor_sites;
        this->acceptor_sites = acceptor_sites;
        this->splice_only = only;
    }
};

#endif // SPLICE_H
