#include <iostream>
#include <fstream>
using namespace std;

int main(int argc, char *argv[]){
    cout<<"hellow"<<endl;
    char hepmc[1000];
    char fastjet[1000];
    cout << "Type hepmc2 absolute directory: "; // Type a number and press enter
    cin >> hepmc; // Get user input from the keyboard
    cout << "Type fastjet3 absolute directory: "; // Type a number and press enter
    cin >> fastjet; // Get user input from the keyboard
    cout << "Your hepmc2 is at: " << hepmc <<"\n"<< "Your fastjet is at: " << fastjet <<endl; // Display the input value
    ofstream myfile;
    myfile.open("config.in", ios::out);
    myfile<<"HEPMC2="<<hepmc<<"\n"<<"FASTJET3="<<fastjet;
    myfile.close();
    return 0;
}