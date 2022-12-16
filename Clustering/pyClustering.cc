#include <vector>
#include "ClusteringAlgorithm.hh"

#include <pybind11.h>
#include <stl.h>
#include <stdint.h>

namespace py = pybind11;

using namespace std;

std::vector<vector<vector<int64_t>>> clustering(int64_t numEvt, std::vector<int64_t> eventsListDirect, std::vector<int64_t> copyListDirect, std::vector<int64_t> strandListDirect, std::vector<int64_t> eventsListIndirect, std::vector<int64_t> copyListIndirect, std::vector<int64_t> strandListIndirect)
{
    std::vector<vector<int64_t>> SBresults;
    std::vector<vector<int64_t>> clusterFrequency; // direct, indirect, total

    // start clustering object
    ClusteringAlgorithm *fpClusteringDirect = new ClusteringAlgorithm(10, 2, 5, 37.5);   // eV
    ClusteringAlgorithm *fpClusteringIndirect = new ClusteringAlgorithm(10, 2, 5, 37.5); // eV
    ClusteringAlgorithm *fpClusteringTotal = new ClusteringAlgorithm(10, 2, 5, 37.5);    // eV

    // Clustering
    for (int64_t e = 0; e < numEvt; ++e)
    {
        std::vector<int64_t> tempResults(13, 0);
        std::vector<int64_t> clusterFrequencyTemp(51, 0); // direct, indirect, total

        bool eventFound =false;

        for (int64_t i = 0; i < eventsListDirect.size(); ++i)
        {
            if (e == eventsListDirect[i])
            {
                fpClusteringDirect->RegisterDamage(copyListDirect[i], strandListDirect[i], 1);
                fpClusteringTotal->RegisterDamage(copyListDirect[i], strandListDirect[i], 1); 
                eventFound=true;
            }
        }

        for (int64_t i = 0; i < eventsListIndirect.size(); ++i)
        {
            if (e == eventsListIndirect[i])
            {
                fpClusteringIndirect->RegisterDamage(copyListIndirect[i], strandListIndirect[i], 2);
                fpClusteringTotal->RegisterDamage(copyListIndirect[i], strandListIndirect[i], 2); 
                eventFound=true;
            }
        }

        if (!eventFound) continue;

        std::map<int64_t, int64_t> sizeDistributionDirect = fpClusteringDirect->RunClustering();
        tempResults[0] = e;
        tempResults[1] = fpClusteringDirect->GetTotalSB(1);
        tempResults[2] = fpClusteringDirect->GetSSB(1);
        tempResults[3] = fpClusteringDirect->GetComplexSSB(1);
        tempResults[4] = fpClusteringDirect->GetDSB(1);


        fpClusteringDirect->Purge();

        sizeDistributionDirect.clear();

        std::map<int64_t, int64_t> sizeDistributionIndirect = fpClusteringIndirect->RunClustering();

        tempResults[5] = fpClusteringIndirect->GetTotalSB(2);
        tempResults[6] = fpClusteringIndirect->GetSSB(2);
        tempResults[7] = fpClusteringIndirect->GetComplexSSB(2);
        tempResults[8] = fpClusteringIndirect->GetDSB(2);


        fpClusteringIndirect->Purge();

        sizeDistributionIndirect.clear();

        std::map<int64_t, int64_t> sizeDistribution = fpClusteringTotal->RunClustering();
        tempResults[9] = fpClusteringTotal->GetTotalSB(0);
        tempResults[10] = fpClusteringTotal->GetSSB(0);
        tempResults[11] = fpClusteringTotal->GetComplexSSB(0);
        tempResults[12] = fpClusteringTotal->GetDSB(0);

        clusterFrequencyTemp[0] = e;
        int c = 1;
        std::map<int64_t, int64_t> sizeDistributionDSB = fpClusteringTotal->GetDSBClusterSizeDistribution(1); //direct

        while ( (!sizeDistributionDSB.empty())&&(sizeDistributionDSB.begin()->first <=10 ))
        {
            clusterFrequencyTemp[c+sizeDistributionDSB.begin()->first-1] = sizeDistributionDSB.begin()->second; //-1 zero indexing
            sizeDistributionDSB.erase(sizeDistributionDSB.begin());
        } 
        
        sizeDistribution.clear();
        sizeDistributionDSB.clear();

        c = 11;
        sizeDistributionDSB = fpClusteringTotal->GetDSBClusterSizeDistribution(2); // indirect

        while ( (!sizeDistributionDSB.empty())&&(sizeDistributionDSB.begin()->first <=10 ))
        {
            clusterFrequencyTemp[c+sizeDistributionDSB.begin()->first-1] = sizeDistributionDSB.begin()->second; //-1 zero indexing
            sizeDistributionDSB.erase(sizeDistributionDSB.begin());
        } 
        
        sizeDistribution.clear();
        sizeDistributionDSB.clear();

        c = 21;
        sizeDistributionDSB = fpClusteringTotal->GetDSBClusterSizeDistribution(3); //hybrid

        while ( (!sizeDistributionDSB.empty())&&(sizeDistributionDSB.begin()->first <=10 ))
        {
            clusterFrequencyTemp[c+sizeDistributionDSB.begin()->first-1] = sizeDistributionDSB.begin()->second; //-1 zero indexing
            sizeDistributionDSB.erase(sizeDistributionDSB.begin());
        } 
        
        sizeDistribution.clear();
        sizeDistributionDSB.clear();

        c = 31;

        sizeDistributionDSB = fpClusteringTotal->GetDSBClusterSizeDistribution(4); //mixed

        while ( (!sizeDistributionDSB.empty())&&(sizeDistributionDSB.begin()->first <=10 ))
        {
            clusterFrequencyTemp[c+sizeDistributionDSB.begin()->first-1] = sizeDistributionDSB.begin()->second; //-1 zero indexing
            sizeDistributionDSB.erase(sizeDistributionDSB.begin());
        } 
        
        sizeDistribution.clear();
        sizeDistributionDSB.clear();

        c = 41;

        sizeDistributionDSB = fpClusteringTotal->GetDSBClusterSizeDistribution(0);

        while ( (!sizeDistributionDSB.empty())&&(sizeDistributionDSB.begin()->first <=10 ))
        {
            clusterFrequencyTemp[c+sizeDistributionDSB.begin()->first-1] = sizeDistributionDSB.begin()->second; //-1 zero indexing
            sizeDistributionDSB.erase(sizeDistributionDSB.begin());

        } 
        
        sizeDistribution.clear();
        sizeDistributionDSB.clear();


        fpClusteringTotal->Purge();
        
        SBresults.push_back(tempResults);
        clusterFrequency.push_back(clusterFrequencyTemp);

    }
    std::vector<vector<vector<int64_t>>> results{SBresults, clusterFrequency};

    delete fpClusteringDirect;
    delete fpClusteringIndirect;
    delete fpClusteringTotal;

    return results;
}

PYBIND11_MODULE(clustering, m)
{
    m.doc() = "clustering"; // optional module docstring

    m.def("clustering", &clustering, "A function that performs clustering");
}