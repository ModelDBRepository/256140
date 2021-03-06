/***************************************************************************
 *                           SinState.cpp                                  *
 *                           -------------------                           *
 * copyright            : (C) 2011 by Jesus Garrido                        *
 * email                : jgarrido@atc.ugr.es                              *
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "../../include/learning_rules/SinState.h"

#include "../../include/simulation/ExponentialTable.h"
#include "../../include/simulation/TrigonometricTable.h"


#include <string.h>

#include <iostream>
#include <cstdio>
using namespace std;

#define A 1.0f/2.0f



const float SinState::terms[11][11]  = {{1,0,0,0,0,0,0,0,0,0,0},
	{A,-A,0,0,0,0,0,0,0,0,0},
	{3./2.*pow(A,2),-4./2.*pow(A,2),1./2.*pow(A,2),0,0,0,0,0,0,0,0},
	{10./4.*pow(A,3),-15./4.*pow(A,3),6./4.*pow(A,3),-1./4.*pow(A,3),0,0,0,0,0,0,0},
	{35./8.*pow(A,4),-56./8.*pow(A,4),28./8.*pow(A,4),-8./8.*pow(A,4),1./8.*pow(A,4),0,0,0,0,0,0},
	{126./16.*pow(A,5),-210./16.*pow(A,5),120./16.*pow(A,5),-45./16.*pow(A,5),10./16.*pow(A,5),-1./16.*pow(A,5),0,0,0,0,0},
	{231./16.*pow(A,6),-99./4.*pow(A,6),495./32.*pow(A,6),-55./8.*pow(A,6),66./32.*pow(A,6),-3./8.*pow(A,6),1./32.*pow(A,6),0,0,0,0},
	{429./16.*pow(A,7),-3003./64.*pow(A,7),1001./32.*pow(A,7),-1001./64.*pow(A,7),91./16.*pow(A,7),-91./64.*pow(A,7),7./32.*pow(A,7),-1./64.*pow(A,7),0,0,0},
	{6435./128.*pow(A,8),-715./8.*pow(A,8),1001./16.*pow(A,8),-273./8.*pow(A,8),455./32.*pow(A,8),-35./8.*pow(A,8),15./16.*pow(A,8),-1./8.*pow(A,8),1./128.*pow(A,8),0,0},
	{12155./128.*pow(A,9),-21879./128.*pow(A,9),1989./16.*pow(A,9),-4641./64.*pow(A,9),1071./32.*pow(A,9),-765./64.*pow(A,9),51./16.*pow(A,9),-153./256.*pow(A,9),9./128.*pow(A,9),-1./256.*pow(A,9),0},
	{46189./256.*pow(A,10),-20995./64.*pow(A,10),62985./256.*pow(A,10),-4845./32.*pow(A,10),4845./64.*pow(A,10),-969./32.*pow(A,10),4845./512.*pow(A,10),-285./128.*pow(A,10),95./256.*pow(A,10),-5./128.*pow(A,10),1./512.*pow(A,10)}};


SinState::SinState(unsigned int NumSynapses, unsigned int NewExponent, float NewMaxpos): ConnectionState(NumSynapses, NewExponent+2), exponent(NewExponent), maxpos(NewMaxpos){

	this->tau = this->maxpos/atan((float)exponent);
	this->factor = 1./(exp(-atan((float)this->exponent))*pow(sin(atan((float)this->exponent)),(int) this->exponent));

	if (this->tau==0){
		this->tau = 1e-6;
	}
	inv_tau=1.0f/tau;

	unsigned int ExponenLine = NewExponent/2;

	TermPointer = this->terms[ExponenLine]; 
}

SinState::~SinState() {
}

float SinState::GetPresynapticActivity(unsigned int index){
	return this->GetStateVariableAt(index, 0);
}

float SinState::GetPostsynapticActivity(unsigned int index){
	return 0.0f;
}

unsigned int SinState::GetNumberOfPrintableValues(){
	return ConnectionState::GetNumberOfPrintableValues()+2;
}

double SinState::GetPrintableValuesAt(unsigned int position){
	if (position<ConnectionState::GetNumberOfPrintableValues()){
		return ConnectionState::GetStateVariableAt(0,position);
	} else if (position==ConnectionState::GetNumberOfPrintableValues()) {
		return this->exponent;
	} else if (position==ConnectionState::GetNumberOfPrintableValues()+1) {
		return this->maxpos;
	} else return -1;
}


void SinState::SetNewUpdateTime (unsigned int index, double NewTime, bool pre_post){
	// Update the activity value
	float OldExpon = this->GetStateVariableAt(index, 1);

	float ElapsedTime = float(NewTime - this->GetLastUpdateTime(index));
	float ElapsedRelative = ElapsedTime*this->inv_tau;

	float expon = ExponentialTable::GetResult(-ElapsedRelative);

	this->SetLastUpdateTime(index, NewTime);

	float NewExpon = OldExpon * expon;
	float NewActivity =NewExpon*TermPointer[0];

	int aux=TrigonometricTable::CalculateOffsetPosition(2*ElapsedRelative);
	int LUTindex=0;

	float SinVar, CosVar, OldVarCos, OldVarSin, NewVarCos, NewVarSin;
	int grade, offset;
	for (grade=2, offset=1; grade<=this->exponent; grade+=2, offset++){

		LUTindex =TrigonometricTable::CalculateValidPosition(LUTindex,aux);

		OldVarCos = this->GetStateVariableAt(index, grade);
		OldVarSin = this->GetStateVariableAt(index, grade + 1);

		SinVar = TrigonometricTable::GetElement(LUTindex);
		CosVar = TrigonometricTable::GetElement(LUTindex+1);

		NewVarCos = (OldVarCos*CosVar-OldVarSin*SinVar)*expon;
		NewVarSin = (OldVarSin*CosVar+OldVarCos*SinVar)*expon;

		NewActivity+= NewVarCos*TermPointer[offset];

		this->SetStateVariableAt(index, grade , NewVarCos, NewVarSin);

	}
	NewActivity*=this->factor;
	this->SetStateVariableAt(index, 0, NewActivity, NewExpon);
}





void SinState::ApplyPresynapticSpike(unsigned int index){
	this->incrementStateVaraibleAt(index, 1, 1.0f);
	for (unsigned int grade=2; grade<=this->exponent; grade+=2){
		this->incrementStateVaraibleAt(index, grade, 1.0f);
	}
}

void SinState::ApplyPostsynapticSpike(unsigned int index){
	return;
}


