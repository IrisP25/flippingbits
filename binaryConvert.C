/* Written by Iris Ponce
The goal is to write a code which can convert decimals to binary, 12-bit.
*/
#include <iostream>
#include <bitset>
#include <vector>

using namespace std;
void binaryConvert()
{
    double bits[12] = {0};
    vector<int> t_bit (12),t_bit2 (12);
    //flipped signal
    int signal = 1004;
    //unflipped signal
    string binary = bitset<12>(499).to_string();
    string binary2 = bitset<12>(signal).to_string();
    string zero = bitset<1>(0).to_string();
    int tare =  (int)zero[0];
    cout << binary << endl;
    cout << binary2 << " " << binary2.length() << endl;
    
    for (int i = 0; i < 12; i++)
    {
        bits[i] = (int)binary[i];
        t_bit[i] = ((int)binary[i]-tare);
        t_bit2[i] = ((int)binary2[i]-tare);
 //       cout << (int)t_bit[i] <<" " << i<< endl;
    }
   // cout << sizeof(bits) << endl;
  //  cout << t_bit.size() << endl;
    for (int i = 0;i < t_bit.size(); i++)
    //for (int i = t_bit.size(); i > -1; i--)
    {
         cout << t_bit[i] << " " << t_bit2[i] <<   endl;
         if ( t_bit[i] != t_bit2[i] && t_bit2[i]==1)
         {
             cout << "flipped: " << 11-i << endl;
         }
        
    }
    //805 545.5 522.5
    //2000 474.5 534.5

    



}
