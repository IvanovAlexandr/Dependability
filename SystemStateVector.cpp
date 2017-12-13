/*
  @file SystemStateVector.cpp
  @brief Implementation computing of reliability
  @author Alexandr Ivanov (ivanov.alexandr.1995@gmail.com)
*/

#include <cmath>
#include <iostream>
#include <algorithm>
#include "SystemStateVector.h"

namespace {
    const uint16_t numberOfD = 6;
    const uint16_t numberOfC = 5;
    const uint16_t numberOfB = 3;
    const uint16_t numberOfM = 1;
    const uint16_t numberOfA = 2;
    const uint16_t numberOfPr = 4;
    const uint16_t numberOfElements = numberOfD + numberOfC
                                      + numberOfB + numberOfM
                                      + numberOfA + numberOfPr;

    const double probabilityPr = 1.2 * pow(10, -4);
    const double probabilityA = 1.1 * pow(10, -4);
    const double probabilityC = 1.9 * pow(10, -4);
    const double probabilityD = 3.2 * pow(10, -5);
    const double probabilityB = 1.4 * pow(10, -5);
    const double probabilityM = 3.3 * pow(10, -4);

    enum Sensors { D1 = numberOfPr + numberOfA + numberOfM + numberOfB + numberOfC,
        D2,
        D3,
        D6,
        D7,
        D8};

    enum Controllers { C1 = numberOfPr + numberOfA + numberOfM + numberOfB,
        C2,
        C4,
        C5,
        C6};

    enum Buses { B1 = numberOfPr + numberOfA + numberOfM,
        B2,
        B3};

    enum Highways { M1 = numberOfPr + numberOfA
    };

    enum Adapters { A1 = numberOfPr,
        A2};

    enum Processors { Pr1 = 0,
        Pr2,
        Pr3,
        Pr6};

    const uint16_t currentLoad[] = {50, 45, 60, 50};
    const uint16_t maxLoad[] = {80, 100, 70, 110};

    std::vector<std::vector<uint16_t>> loadRedistributionTab;
    std::vector<uint32_t > failedElements;

    enum State {
    	PASSIVE,
    	ACTIVE
    };

    void initFailedElements() {
        for (uint8_t i = 0; i < numberOfElements; ++i) {
            failedElements.push_back(0);
        }
    }

    void printFailedElements() {
        std::cout << "Failed Elements : "<< std::endl;
        for (int i = 0; i < numberOfElements; ++i) {
            std::cout << "\t" << i <<" : \t" << failedElements[i] << std::endl;
        }
    }

    void createLoadRedistributionTab() {
        for (uint16_t row = 0; row < numberOfPr; ++row) {
            loadRedistributionTab.emplace_back();
            for (uint16_t col = 0; col < numberOfPr; ++col) {
                if (row != col) {
                    loadRedistributionTab.back().push_back((currentLoad[row] < maxLoad[col] - currentLoad[col])
                                                           ? currentLoad[row] : maxLoad[col] - currentLoad[col]);
                } else {
                    loadRedistributionTab.back().push_back(0);
                }
            }
        }
    }

    void printLoadRedistributionTab() {
        for (uint16_t i = 0; i < numberOfPr; ++i) {
            for (uint16_t j = 0; j < numberOfPr; ++j) {
                std::cout << loadRedistributionTab[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    void updateFailedElements(bool f1, bool f2, bool f3, bool f4, std::vector<uint8_t> &vccVec) {
        if (!f1) {
            failedElements[D1] += vccVec[D1];
            failedElements[D2] += vccVec[D2];
            failedElements[C1] += vccVec[C1];
            failedElements[B1] += vccVec[B1];
            failedElements[B2] += vccVec[B2];
            failedElements[Pr1] += vccVec[Pr1];
            failedElements[A1] += vccVec[A1];
            failedElements[M1] += vccVec[M1];
            //failedElements[M2] += vccVec[M2];
            failedElements[A2] += vccVec[A2];
            failedElements[B3] += vccVec[B3];
            //failedElements[B4] += vccVec[B4];
            failedElements[Pr3] += vccVec[Pr3];
        }
        if (!f2) {
            failedElements[D2] += vccVec[D2];
            failedElements[D3] += vccVec[D3];
            failedElements[C2] += vccVec[C2];
            failedElements[B1] += vccVec[B1];
            failedElements[B2] += vccVec[B2];
            failedElements[Pr2] += vccVec[Pr2];
            failedElements[A1] += vccVec[A1];
            failedElements[M1] += vccVec[M1];
            //failedElements[M2] += vccVec[M2];
            failedElements[A2] += vccVec[A2];
            failedElements[B3] += vccVec[B3];
            //failedElements[B4] += vccVec[B4];
            failedElements[Pr6] += vccVec[Pr6];
        }
        if (!f3) {
            failedElements[D7] += vccVec[D7];
            failedElements[C5] += vccVec[C5];
            failedElements[D8] += vccVec[D8];
            failedElements[C6] += vccVec[C6];
            failedElements[B3] += vccVec[B3];
            //failedElements[B4] += vccVec[B4];
            failedElements[Pr6] += vccVec[Pr6];
        }
        if (!f4) {
            failedElements[D6] += vccVec[D6];
            failedElements[C4] += vccVec[C4];
            failedElements[M1] += vccVec[M1];
            //failedElements[M2] += vccVec[M2];
            failedElements[A2] += vccVec[A2];
            failedElements[C5] += vccVec[C5];
            failedElements[B3] += vccVec[B3];
            //failedElements[B4] += vccVec[B4];
            failedElements[Pr3] += vccVec[Pr3];
        }
    }

    bool logicalStructureFunction(uint32_t vcc) {
        std::vector<uint8_t> vccVec;
        for (uint32_t i = 0; i < pow(2, numberOfElements); ++i) {
            vccVec.push_back(static_cast<uint8_t>(vcc & 0x1));
            vcc >>= 1;
        }
        /*

        bool f1 = vccVec[D1] && vccVec[D2] && vccVec[C1] && (vccVec[B1] || vccVec[B2])
                && (vccVec[Pr1] || (vccVec[A1] && (vccVec[M1] || vccVec[M2]) && vccVec[A2]
                && (vccVec[B3] || vccVec[B4]) && vccVec[Pr3]));

         bool f2 = vccVec[D2] && vccVec[D3] && vccVec[C2] && (vccVec[B1] || vccVec[B2])
                && (vccVec[Pr2] || (vccVec[A1] && (vccVec[M1] || vccVec[M2]) && vccVec[A2]
                && (vccVec[B3] || vccVec[B4]) && vccVec[Pr6]));

        bool f3 = ((vccVec[D7] && vccVec[C5]) || (vccVec[D8] && vccVec[C6])) && (vccVec[B3] || vccVec[B4]) && vccVec[Pr6];

        bool f4 = vccVec[D6] && (vccVec[C4] && (vccVec[M1] || vccVec[M2]) && vccVec[A2] || vccVec[C5])
                  && (vccVec[B3] || vccVec[B4]) && vccVec[Pr3];

        */
        bool f1 = vccVec[D1] && vccVec[D2] && vccVec[C1] && (vccVec[B1] || vccVec[B2])
                  && (vccVec[Pr1] || (vccVec[A1] && vccVec[M1] && vccVec[A2] && vccVec[B3] && vccVec[Pr3]));
        bool f2 = vccVec[D2] && vccVec[D3] && vccVec[C2] && (vccVec[B1] || vccVec[B2])
                  && (vccVec[Pr2] || (vccVec[A1] && vccVec[M1] && vccVec[A2] && vccVec[B3] && vccVec[Pr6]));
        bool f3 = ((vccVec[D7] && vccVec[C5]) || (vccVec[D8] && vccVec[C6])) && vccVec[B3] && vccVec[Pr6];
        bool f4 = vccVec[D6] && (vccVec[C4] && vccVec[M1] && vccVec[A2] || vccVec[C5]) && vccVec[B3] && vccVec[Pr3];



        updateFailedElements(f1, f2, f3, f4, vccVec);

        return f1 && f2 && f3 && f4;
    }

    double calcProbability(std::vector<uint8_t> vccVec) {
        double probability = 1;

        uint16_t loopStart = 0;
        for (uint16_t i = loopStart; i < numberOfPr; ++i) {
            probability *= vccVec[i] * (1 - probabilityPr) + (1 - vccVec[i]) * probabilityPr;
        }
        loopStart += numberOfPr;
        for (uint16_t i = loopStart; i < loopStart + numberOfA; ++i) {
            probability *= vccVec[i] * (1 - probabilityA) + (1 - vccVec[i]) * probabilityA;
        }
        loopStart += numberOfA;
        for (uint16_t i = loopStart; i < loopStart + numberOfM; ++i) {
            probability *= vccVec[i] * (1 - probabilityM) + (1 - vccVec[i]) * probabilityM;
        }
        loopStart += numberOfM;
        for (uint16_t i = loopStart; i < loopStart + numberOfB; ++i) {
            probability *= vccVec[i] * (1 - probabilityB) + (1 - vccVec[i]) * probabilityB;
        }
        loopStart += numberOfB;
        for (uint16_t i = loopStart; i < loopStart + numberOfC; ++i) {
            probability *= vccVec[i] * (1 - probabilityC) + (1 - vccVec[i]) * probabilityC;
        }
        loopStart += numberOfC;
        for (uint16_t i = loopStart; i < loopStart + numberOfD; ++i) {
            probability *= vccVec[i] * (1 - probabilityD) + (1 - vccVec[i]) * probabilityD;
        }

        return probability;
    }

    double calcFaultTolerance(std::vector<std::pair<uint32_t, uint32_t>> &vcc, double size, State state) {
        double faultTolerance = 0;
    	for (uint16_t vccId = 0; vccId < vcc.size() * size; ++vccId) {
    		uint32_t targetVcc = (state == PASSIVE) ? vcc[vccId].first : vcc[vccId].second;

	        if (logicalStructureFunction(targetVcc)) {
	            uint32_t tempVcc = vcc[vccId].first;
	            std::vector<uint8_t> vccVec;
	            for (uint32_t i = 0; i < pow(2, numberOfElements); ++i) {
	                vccVec.push_back(static_cast<uint8_t>(tempVcc & 0x1));
	                tempVcc >>= 1;
	            }

	            faultTolerance += calcProbability(vccVec);
	        }	
        
    	}
        return faultTolerance;
    }
}


void SystemStateVector::generateVcc() {
    for (uint32_t i = 1; i < pow(2, numberOfElements); ++i) {
        uint32_t vcc = i;
        uint8_t numberOfZero = 0;
        for (uint8_t bit = 0; bit < numberOfElements; ++bit) {
            numberOfZero += ((vcc & 0x1) == 0) ? 1 : 0;
            if (numberOfZero > 4) {
                break;
            }
            vcc >>= 1;
        }
        switch (numberOfZero) {
            case 1: {
                m_vcc1zero.emplace_back(i, i);
                break;
            }
            case 2: {
                m_vcc2zero.emplace_back(i, i);
                break;
            }
            case 3: {
                m_vcc3zero.emplace_back(i, i);
                break;
            }
            case 4: {
                m_vcc4zero.emplace_back(i, i);
                break;
            }
            default:break;
        }
    }


}

void SystemStateVector::generateVirtualSystemStateVector() {

    createLoadRedistributionTab();
    printLoadRedistributionTab();
    createVirtualSystemStateVector(m_vcc1zero);
    createVirtualSystemStateVector(m_vcc2zero);
    createVirtualSystemStateVector(m_vcc3zero);
    createVirtualSystemStateVector(m_vcc4zero);
}

void SystemStateVector::createVirtualSystemStateVector(std::vector<std::pair<uint32_t, uint32_t>> &realVcc) {
    for (auto &item : realVcc) {
        uint32_t vcc = item.first;
        std::vector<uint8_t> processors;
        for (uint8_t bit = 0; bit < loadRedistributionTab.size(); ++bit) {
            processors.push_back(static_cast<uint8_t>(vcc & 0x1));
            vcc >>= 1;
        }

        std::vector<std::vector<uint16_t>> currentLoadRedistributionTab = loadRedistributionTab;
        for (uint16_t i = 0; i < processors.size(); ++i) {
            if (processors[i] == 0) {
                uint16_t load = currentLoad[i];
                for (uint16_t j = 0; j < processors.size(); ++j) {
                    if (j == i || processors[j] == 0) continue;
                    if (currentLoadRedistributionTab[i][j] > load) {
                        for (auto& row : currentLoadRedistributionTab) {
                            row[j] = (row[j] > load) ? row[j] - load : load - row[j];
                        }
                        load = 0;
                    } else {
                        uint16_t redistribution = currentLoadRedistributionTab[i][j];
                        for (auto& row : currentLoadRedistributionTab) {
                            row[j] = (row[j] > load) ? row[j] - redistribution : row[j] - row[j];
                        }
                        load -= redistribution;
                    }
                    if (load == 0) break;
                }
                if (load == 0) item.second = item.second | 1 << i;

            }
        }
    }

}

std::vector<double> SystemStateVector::getPassiveFaultTolerance() {
    initFailedElements();
    std::vector<double> allFaultTolerance;
    std::vector<uint8_t> vec;
    for (uint8_t i = 0; i < numberOfElements; ++i) {
        vec.push_back(1);
    }

    allFaultTolerance.push_back(calcProbability(vec));

    allFaultTolerance.push_back(calcFaultTolerance(m_vcc1zero, 1, PASSIVE));

    allFaultTolerance.push_back(calcFaultTolerance(m_vcc2zero, 1, PASSIVE));

    std::shuffle( m_vcc3zero.begin(), m_vcc3zero.end(), std::mt19937(std::random_device()()));
    allFaultTolerance.push_back(calcFaultTolerance(m_vcc3zero, 0.5, PASSIVE) * 2);

    std::shuffle( m_vcc4zero.begin(), m_vcc4zero.end(), std::mt19937(std::random_device()()));
    allFaultTolerance.push_back(calcFaultTolerance(m_vcc4zero, 0.1, PASSIVE) * 10);

    return allFaultTolerance;
}

std::vector<double> SystemStateVector::getActiveFaultTolerance() {

    initFailedElements();
    std::vector<double> allFaultTolerance;
    std::vector<uint8_t> vec;
    for (uint8_t i = 0; i < numberOfElements; ++i) {
        vec.push_back(1);
    }

    allFaultTolerance.push_back(calcProbability(vec));

    allFaultTolerance.push_back(calcFaultTolerance(m_vcc1zero, 1, ACTIVE));

    allFaultTolerance.push_back(calcFaultTolerance(m_vcc2zero, 1, ACTIVE));

    std::shuffle( m_vcc3zero.begin(), m_vcc3zero.end(), std::mt19937(std::random_device()()));
    allFaultTolerance.push_back(calcFaultTolerance(m_vcc3zero, 0.5, ACTIVE) * 2);

    std::shuffle( m_vcc4zero.begin(), m_vcc4zero.end(), std::mt19937(std::random_device()()));
    allFaultTolerance.push_back(calcFaultTolerance(m_vcc4zero, 0.1, ACTIVE) * 10);

    printFailedElements();

    return allFaultTolerance;
}


unsigned long SystemStateVector::getSizeVcc1zero() {
    return m_vcc1zero.size();
}

unsigned long SystemStateVector::getSizeVcc2zero() {
    return m_vcc2zero.size();
}

unsigned long SystemStateVector::getSizeVcc3zero() {
    return m_vcc3zero.size();
}

unsigned long SystemStateVector::getSizeVcc4zero() {
    return m_vcc4zero.size();
}
