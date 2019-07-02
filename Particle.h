//  cHNLdecay  --  calculate decay widths of Heavy Neutral Leptons
//  Copyright (C) 2018 - Fabian A.J. Thiele, <fabian.thiele@posteo.de>
//
//  This file is part of cHNLdecay.
//
//  cHNLdecay is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  cHNLdecay is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef PARTICLE_H
#define PARTICLE_H

#include "TString.h"
#include <map>

class Particle {
public:
  Particle() {
    pdgid = 0;
    name = pdgIdToLaTeX(0);
    mass = 0;
  }
  Particle(Int_t p, Double_t m) {
    pdgid = p;
    name = pdgIdToLaTeX(p);
    mass = m;
  }
  Particle(const Particle &obj) {
    pdgid = obj.getPdgId();
    name = obj.getName();
    mass = obj.getMass();
  }

  Double_t getMass() const { return mass; }
  TString getName() const { return name; }

  Int_t getPdgId() const { return pdgid; }

  bool operator==(const Particle &a) const {
    return pdgid == a.getPdgId(); // true if they have the same pdg IDs
  }

  Particle &operator=(const Particle &obj) {
    pdgid = obj.getPdgId();
    name = obj.getName();
    mass = obj.getMass();

    return *this;
  }

protected:
  TString pdgIdToLaTeX(Int_t p) const {
    std::map<Int_t, TString> label = {
        {1, "d"}, {2, "u"}, {3, "s"}, {4, "c"},   {5, "b"}, {6, "t"},
        {11, "e^-"},           	{12, "\\nu_e"},   {13, "\\mu^-"},
        {14, "\\nu_\\mu"},     	{15, "\\tau^-"},  {16, "\\nu_\\tau"},
        {111, "\\pi^0"},       	{211, "\\pi^0"},  {113, "\\rho^0"},   {213, "\\rho^+"}, 
        {311, "K^0"},			{321, "K^+"},
        {313, "K^{\\ast0}"},	{323, "K^{\\ast+}"},
        {411, "D^+"},			{421, "D^0"},	  	 {431, "D_s^+"}, 
        {413, "D^{\\ast+}"},	{423, "D^{\\ast+}"}, {433, "D^{\\ast+}_s"},
        {511, "B^0"},			{521, "B^+"},     	 {531, "B_s^0"},		{541, "B_c^+"},
		{221, "\\eta"},   		{331, "\\eta'"},
        {441, "\\eta_c"},      	{213, "\\rho^+"}, 
		{223, "\\omega"},		{333, "\\phi"},      {443, "J/\\Psi"}};

    if (p > 0)
      return label.at(p);
    else if (p < 0) {
      TString ret = label.at(-p);
      if (ret.Contains("+"))
        ret.ReplaceAll("+", "-");
      else if (ret.Contains("-"))
        ret.ReplaceAll("-", "+");
      else
        ret = TString("\\bar{" + ret + "}");

      return ret;
    }
    return TString("");
  }

  Int_t pdgid;
  TString name;
  Double_t mass;
};
#endif
