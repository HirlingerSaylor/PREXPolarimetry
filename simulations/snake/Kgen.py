#!/usr/bin/python

from ROOT import TCanvas, TH1F, TF1
import math, random

def Kgen():

    kh1   = TH1F("kh1", "kh1", 200, 0., 0.0002)
    kh2   = TH1F("kh2", "kh2", 200, 0., 0.0002)
    mymax = 0.00526615620525
    me    = 0.000511
    c10 = 1.0
    Z10 = 26.0 - 0.0 - ( 2.0 - 1.0 ) / 2.0
    P10 = Z10
    #n   = 10000000
    #for i in range(n):
    while 1 == 1:
        har   = random.random() * 500.
        p     = har
        p_nat = har * me / 137.
        prob  = random.random() * mymax
        
        q   = p / P10
        ans = 4.0 * c10 * q ** 2 / math.pi / P10 / ( 1.0 + q ** 2.0 ) ** 4.0
        
        if prob <= ans:
            #print  p_nat
            return p_nat
        #kh1.Fill(p_nat)
        #kh2.Fill(p_nat, ans)
    #can = TCanvas("can", "can", 20, 20, 900, 900)
    #can.SetLogy()
    #kh1.Scale(1. / kh1.Integral(0, 200) )
    #kh1.GetYaxis().SetRangeUser(0.00001, 1.)
    #kh2.Scale(1. / kh2.Integral(0, 200) )
    #kh2.GetYaxis().SetRangeUser(0.00001, 1.)


    #f1  = TF1("f1", "4.0 * x * x / 0.000511 / 0.000511 * 137. * 137. / 25.5 / 25.5 / 3.141592654 / 25.5 / TMath::Power( ( 1.0 + x * x / 0.000511 / 0.000511 * 137. * 137. / 25.5 / 25.5 ) , 4.0 )", 0., 0.002)
    #norm = f1.Integral(0., 0.002) / 0.002 * 200.
    #print norm
    #f2  = TF1("f2", "4.0 * x * x / 0.000511 / 0.000511 * 137. * 137. / 25.5 / 25.5 / 3.141592654 / 25.5 / TMath::Power( ( 1.0 + x * x / 0.000511 / 0.000511 * 137. * 137. / 25.5 / 25.5 ) , 4.0 ) / " + str(norm), 0., 0.002)
    #f2.Draw()
    #kh1.Draw()
    #kh2.Draw("SAME")
    press = raw_input("Press enter to continue.")

if __name__ == '__main__':
    Kgen()
