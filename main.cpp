#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
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



}



int main() {

    std::srand ( unsigned ( std::time(0) ) );


    SystemStateVector systemStateVector;
    systemStateVector.generateVcc();
    systemStateVector.generateVirtualSystemStateVector();
    //createLoadRedistributionTab();


    std::cout << "1 zero: " << systemStateVector.getSizeVcc1zero() << std::endl;
    //printTab(vcc1zero);
    std::cout << "2 zero: " << systemStateVector.getSizeVcc2zero() << std::endl;
    //printTab(vcc2zero);
    std::cout << "3 zero: " << systemStateVector.getSizeVcc3zero() << std::endl;
    //printTab(vcc3zero);
    std::cout << "4 zero: " << systemStateVector.getSizeVcc4zero() << std::endl;
    //printTab(vcc4zero);



    std::cout << "Passive Fault Tolerance = " << std::endl;
    std::vector<double> passiveFaultTolerance = systemStateVector.getPassiveFaultTolerance();

    std::cout << "0 fails :\t" << passiveFaultTolerance[0] << std::endl;
    std::cout << "1 fails :\t" << passiveFaultTolerance[1] << std::endl;
    std::cout << "2 fails :\t" << passiveFaultTolerance[2] << std::endl;
    std::cout << "3 fails :\t" << passiveFaultTolerance[3] << std::endl;
    std::cout << "4 fails :\t" << passiveFaultTolerance[4] << std::endl;

    std::cout << "Work probability = \t" << passiveFaultTolerance[0] +
                                            passiveFaultTolerance[1] +
                                            passiveFaultTolerance[2] +
                                            passiveFaultTolerance[3] +
                                            passiveFaultTolerance[4] << std::endl;

    std::cout << "Fail probability = \t" << 1 - (passiveFaultTolerance[0] +
                                            passiveFaultTolerance[1] +
                                            passiveFaultTolerance[2] +
                                            passiveFaultTolerance[3] +
                                            passiveFaultTolerance[4]) << std::endl;

   // for (uint8_t i = 0; i < numberOfElements; ++i) {
   //     failedElements[i] = 0;
   // }

    std::cout << "\nActive Fault Tolerance = " << std::endl;
    std::vector<double> activeFaultTolerance = systemStateVector.getActiveFaultTolerance();

    std::cout << "0 fails :\t" << activeFaultTolerance[0] << std::endl;
    std::cout << "1 fails :\t" << activeFaultTolerance[1] << std::endl;
    std::cout << "2 fails :\t" << activeFaultTolerance[2] << std::endl;
    std::cout << "3 fails :\t" << activeFaultTolerance[3] << std::endl;
    std::cout << "4 fails :\t" << activeFaultTolerance[4] << std::endl;

    std::cout << "Work probability = \t" << activeFaultTolerance[0] +
                                            activeFaultTolerance[1] +
                                            activeFaultTolerance[2] +
                                            activeFaultTolerance[3] +
                                            activeFaultTolerance[4] << std::endl;

    std::cout << "Fail probability = \t" << 1 - (activeFaultTolerance[0] +
                                                 activeFaultTolerance[1] +
                                                 activeFaultTolerance[2] +
                                                 activeFaultTolerance[3] +
                                                 activeFaultTolerance[4]) << std::endl;


   // std::cout << "Failed Elements : "<< std::endl;
 //   for (int i = 0; i < numberOfElements; ++i) {
 //       std::cout << "\t" << i <<" : \t" << failedElements[i]*20 << std::endl;
 //   }

    //std::cout << "LoadRedistributionTab: " << std::endl;
    //printLoadRedistributionTab();
    std::cout << "Press enter to continue ...";
    //std::cin.get();

    return 0;
}