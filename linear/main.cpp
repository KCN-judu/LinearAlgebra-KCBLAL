#include "Linear.h"
#include <iostream>
using namespace std;
int main(){
    pair<int,int> half = make_pair(1,2);
    cout<<KLinear::add(half,half).first<<"/"<<KLinear::add(half,half).second;
    getchar();
}