#ifndef CALLBACK_H
#define CALLBACK_H

#include <vector>

class NomadLinker;

class Callback
{
    public:
        Callback() {}
        virtual ~Callback() {}
        virtual void call( NomadLinker& object, std::vector<double> xvec ) {}
};

#endif
