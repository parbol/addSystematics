//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                              //
//                                                                                              //
// To compile type: make                                                                        //
//                                                                                              //
// To execute type: ./AddSystematics --input sourcefile --output destinyfile --formula file.txt //
//                                                                                              //
//                                                                                              //
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>
#include <fstream>
#include <TFormula.h>

//-----------------------------------------------------------------------//
//Variables of the tree                                                  //
//-----------------------------------------------------------------------//
Float_t         philep1, etalep1, ptlep1, mlep1;
Int_t           idlep1;
Float_t         philep2, etalep2, ptlep2, mlep2;
Int_t           idlep2;
Float_t         phib1, etab1, ptb1, mb1;
Float_t         phib2, etab2, ptb2, mb2;
Float_t         phinu1, etanu1, ptnu1, mnu1;
Float_t         phinu2, etanu2, ptnu2, mnu2;
Float_t         phiChi1, etaChi1, ptChi1, mChi1;
Float_t         phiMET, etaMET, ptMET, mMET;
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//


//-----------------------------------------------------------------------//
//Function declaratione                                                  //
//-----------------------------------------------------------------------//
bool getOptions(int , char **, std::string &, std::string &, std::string &);

void showBanner();

void setVariables(TTree *);

void printError(std::string);

void printLog(std::string);

void applyLeptonSystematic(TRandom3 *, TFormula *, TFormula *);

void applybjetSystematic(TRandom3 *, TFormula *, TFormula *);

void applyMETSystematic(TRandom3 *, TFormula*, TFormula *);
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//


//-----------------------------------------------------------------------//
//Main program                                                           //
//-----------------------------------------------------------------------//
int main(int argc, char **argv) {

    showBanner();


    //------------------------------------------
    //Getting options
    std::string inputFile;
    std::string outputFile;
    std::string formulaFile;
    if(!getOptions(argc, argv, inputFile, outputFile, formulaFile)) {
        printError("Usage: ./addSystematics --input inputFile.root --output outputFile.root --formula formulaFile.txt");
        return -1;
    }
    
    //------------------------------------------
    //Getting input file and tree    
    TFile *fin = new TFile(inputFile.c_str());
    if(fin == NULL) {
        printError("The source file does not exist");
        return -1;
    }
    TTree *tin = (TTree *) fin->Get("t");
    setVariables(tin);

    printLog(std::string("File ") + inputFile + std::string(" successfully opened"));

    
    //------------------------------------------
    //Getting output file and tree
    TFile *fout = new TFile(outputFile.c_str(), "RECREATE");
    if(fout == NULL) {
        printError("The output file cannot be created");
        return -1;
    }
    fout->cd();
    TTree *tout = (TTree *) tin->Clone("t");
    tout->Reset();
    setVariables(tout);
    
    printLog(std::string("File ") + outputFile + std::string(" successfully opened"));


    //------------------------------------------
    //Getting the formulas for the systematics
    std::vector<std::string> formulaNames;
    formulaNames.push_back("Lepton mu");
    formulaNames.push_back("Lepton sigma");
    formulaNames.push_back("b-jet mu");
    formulaNames.push_back("b-jet sigma");
    formulaNames.push_back("met mu");
    formulaNames.push_back("met sigma");
    std::vector<std::string> vector_formula;
    std::vector<TFormula> formulas;
    std::ifstream forma(formulaFile.c_str(), std::ifstream::in);
    int counter = 0;
    while(!forma.eof()) {
        std::string line;
        getline(forma, line);
        if(line.c_str()[0] != 'f') continue;
        if(counter == 6) {
            printError("There was a problem with the formula file");
            return -1;
        }
        std::string formula = line.substr(2);
        vector_formula.push_back(formula);
        TFormula fo(formulaNames[counter].c_str(), formula.c_str());
        if(fo.Compile() != 0) {
            printError(std::string("There was a problem with formula: ") + formula);
            return -1;
        }
        formulas.push_back(fo);
        counter++;
    }
    forma.close();
  
    if(vector_formula.size() != 6) {
        printError("There was a problem with the formula file");
        return -1;
    }

    printLog("Using the following formulas:");
    for(unsigned int i = 0; i < vector_formula.size(); i++) {
        printLog(formulaNames[i] + std::string(":") + vector_formula[i]);
    }

    //------------------------------------------
    //Random number generator
    TRandom3 *ran = new TRandom3();
 

    //------------------------------------------
    //Looping over the tree
    Long64_t nentries = tin->GetEntriesFast();
    for(Long64_t jentry=0; jentry < nentries; jentry++) {
      //Reading entry for i 
      tin->GetEntry(jentry);
      
      //Systematics
      applyLeptonSystematic(ran, &formulas[0], &formulas[1]);

      applybjetSystematic(ran, &formulas[2], &formulas[3]);

      applyMETSystematic(ran, &formulas[4], &formulas[5]);

      //Filling output tree
      tout->Fill();
   }

   tout->Write();
   fout->Close();
   fin->Close();


}
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//


//-----------------------------------------------------------------------//
// Apply systematics to the leptons and MET                              //
//-----------------------------------------------------------------------//
void applyLeptonSystematic(TRandom3 *ran, TFormula *mu, TFormula *sigma) {

    //Evaluating the formulas
    Float_t mu1, sigma1;
    Float_t mu2, sigma2;

    mu1 = mu->Eval(ptlep1, etalep1);
    mu2 = mu->Eval(ptlep2, etalep2);
    sigma1 = sigma->Eval(ptlep1, etalep1);
    sigma2 = sigma->Eval(ptlep2, etalep2);

    Float_t pt1_add = ran->Gaus(mu1, sigma1);
    Float_t pt2_add = ran->Gaus(mu2, sigma2);
    
    //updating leptons
    Float_t pt1_x = (pt1_add - ptlep1) * cos(philep1);
    Float_t pt1_y = (pt1_add - ptlep1) * sin(philep1);
    Float_t pt2_x = (pt2_add - ptlep2) * cos(philep2);
    Float_t pt2_y = (pt2_add - ptlep2) * sin(philep2);

    ptlep1 = pt1_add;
    ptlep2 = pt2_add;

    //Updating MET
    Float_t met_x = (ptMET * cos(phiMET)) - pt1_x - pt2_x; 
    Float_t met_y = (ptMET * sin(phiMET)) - pt1_y - pt2_y;
    ptMET = sqrt(met_x*met_x + met_y*met_y);
    phiMET = atan2(met_y, met_x);


}
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//


//-----------------------------------------------------------------------//
// Apply systematics to the bjets and MET                                //
//-----------------------------------------------------------------------//
void applybjetSystematic(TRandom3 *ran, TFormula *mu, TFormula *sigma) {

    //Evaluating the formulas
    Float_t mu1, sigma1;
    Float_t mu2, sigma2;

    mu1 = mu->Eval(ptb1, etab1);
    mu2 = mu->Eval(ptb2, etab2);
    sigma1 = sigma->Eval(ptb1, etab1);
    sigma2 = sigma->Eval(ptb2, etab2);

    Float_t pt1_add = ran->Gaus(mu1, sigma1);
    Float_t pt2_add = ran->Gaus(mu2, sigma2);
    
    //updating leptons
    Float_t pt1_x = (pt1_add - ptb1) * cos(phib1);
    Float_t pt1_y = (pt1_add - ptb1) * sin(phib1);
    Float_t pt2_x = (pt2_add - ptb2) * cos(phib2);
    Float_t pt2_y = (pt2_add - ptb2) * sin(phib2);

    ptb1 = pt1_add;
    ptb2 = pt2_add;

    //Updating MET
    Float_t met_x = (ptMET * cos(phiMET)) - pt1_x - pt2_x; 
    Float_t met_y = (ptMET * sin(phiMET)) - pt1_y - pt2_y;
    ptMET = sqrt(met_x*met_x + met_y*met_y);
    phiMET = atan2(met_y, met_x);

}
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//


//-----------------------------------------------------------------------//
// Apply systematics to the MET (Unclustered)                            //
//-----------------------------------------------------------------------//
void applyMETSystematic(TRandom3 *ran, TFormula *mu, TFormula *sigma) {

    //Evaluating the formulas
    Float_t mum, sigmam;

    mum = mu->Eval(ptMET);
    sigmam = sigma->Eval(ptMET);

    Float_t met_x = ran->Gaus(mum, sigmam);
    Float_t met_y = ran->Gaus(mum, sigmam);
    
    //Updating MET
    ptMET = sqrt(met_x*met_x + met_y*met_y);
    phiMET = atan2(met_y, met_x);

}
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//


//-----------------------------------------------------------------------//
//Allocate the tree branch into local variables                          //
//-----------------------------------------------------------------------//
void setVariables(TTree *tree) {

   tree->SetBranchAddress("philep1", &philep1);
   tree->SetBranchAddress("etalep1", &etalep1);
   tree->SetBranchAddress("ptlep1", &ptlep1);
   tree->SetBranchAddress("mlep1", &mlep1);
   tree->SetBranchAddress("idlep1", &idlep1);
   tree->SetBranchAddress("philep2", &philep2);
   tree->SetBranchAddress("etalep2", &etalep2);
   tree->SetBranchAddress("ptlep2", &ptlep2);
   tree->SetBranchAddress("mlep2", &mlep2);
   tree->SetBranchAddress("idlep2", &idlep2);
   tree->SetBranchAddress("phib1", &phib1);
   tree->SetBranchAddress("etab1", &etab1);
   tree->SetBranchAddress("ptb1", &ptb1);
   tree->SetBranchAddress("mb1", &mb1);
   tree->SetBranchAddress("phib2", &phib2);
   tree->SetBranchAddress("etab2", &etab2);
   tree->SetBranchAddress("ptb2", &ptb2);
   tree->SetBranchAddress("mb2", &mb2);
   tree->SetBranchAddress("phinu1", &phinu1);
   tree->SetBranchAddress("etanu1", &etanu1);
   tree->SetBranchAddress("ptnu1", &ptnu1);
   tree->SetBranchAddress("mnu1", &mnu1);
   tree->SetBranchAddress("phinu2", &phinu2);
   tree->SetBranchAddress("etanu2", &etanu2);
   tree->SetBranchAddress("ptnu2", &ptnu2);
   tree->SetBranchAddress("mnu2", &mnu2);
   tree->SetBranchAddress("phiChi1", &phiChi1);
   tree->SetBranchAddress("etaChi1", &etaChi1);
   tree->SetBranchAddress("ptChi1", &ptChi1);
   tree->SetBranchAddress("mChi1", &mChi1);
   tree->SetBranchAddress("phiMET", &phiMET);
   tree->SetBranchAddress("etaMET", &etaMET);
   tree->SetBranchAddress("ptMET", &ptMET);
   tree->SetBranchAddress("mMET", &mMET);

}
//-----------------------------------------------------------------------//
//-----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// This method will put the command line input in the variables.        //
//----------------------------------------------------------------------//
bool getOptions(int argc, char **argv, std::string &nameOfInputFile, std::string &nameOfOutputFile, std::string &formulaFile) {

    int option_iterator;
    int option_counter = 0;
    bool moreoptions = true;

    while (moreoptions) {
        static struct option long_options[] = {
            /* These options set a flag. */
            {"input",     required_argument, 0, 'i'},
            {"output",    required_argument, 0, 'o'},
            {"formula",   required_argument, 0, 'f'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        option_iterator = getopt_long(argc, argv, "d:", long_options, &option_index);
        if (option_iterator == -1) {
            moreoptions = false;
        } else {
            option_counter++;
            switch (option_iterator) {
            case 0:
                if (long_options[option_index].flag != 0)
                    break;
                if (optarg)
                    break;
            case 'i':
                nameOfInputFile = (std::string) optarg;
                break;
            case 'o':
                nameOfOutputFile = (std::string) optarg;
                break;
            case 'f':
                formulaFile = (std::string) optarg;
                break;
            case '?':
                return false;
                break;
            default:
                return false;
            }
        }
    }

    if (option_counter == 0) {
        return false;
    }
    return true;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Print errors                                                         //
//----------------------------------------------------------------------//
void printError(std::string error)  {
    std::cout << "\033[1;31m" << "[Error] " << error << "\033[0m" << std::endl << std::endl;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Print Logs                                                         //
//----------------------------------------------------------------------//
void printLog(std::string error)  {
    std::cout << "\033[1;33m" << "[Log] " << "\033[0m" << error << std::endl << std::endl;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//


//----------------------------------------------------------------------//
// Banner                                                               //
//----------------------------------------------------------------------//
void showBanner()  {
     std::cout << std::endl << std::endl;
     std::cout << "\033[1;34m"
               << "     ___       _______   _______       _______.____    ____  _______.___________. _______ .___  ___.      ___   .___________. __    ______     _______." << std::endl
               << "    /   \\     |       \\ |       \\     /       |\\   \\  /   / /       |           ||   ____||   \\/   |     /   \\  |           ||  |  /      |   /       |" << std::endl
               << "   /  ^  \\    |  .--.  ||  .--.  |   |   (----` \\   \\/   / |   (----`---|  |----`|  |__   |  \\  /  |    /  ^  \\ `---|  |----`|  | |  ,----'  |   (----`" << std::endl
               << "  /  /_\\  \\   |  |  |  ||  |  |  |    \\   \\      \\_    _/   \\   \\       |  |     |   __|  |  |\\/|  |   /  /_\\  \\    |  |     |  | |  |        \\   \\    " << std::endl
               << " /  _____  \\  |  '--'  ||  '--'  |.----)   |       |  | .----)   |      |  |     |  |____ |  |  |  |  /  _____  \\   |  |     |  | |  `----.----)   |   " << std::endl
               << "/__/     \\__\\ |_______/ |_______/ |_______/        |__| |_______/       |__|     |_______||__|  |__| /__/     \\__\\  |__|     |__|  \\______|_______/    " << std::endl << "\033[0m" << std::endl;   
     std::cout << std::endl << std::endl;
}
//----------------------------------------------------------------------//
//----------------------------------------------------------------------//

