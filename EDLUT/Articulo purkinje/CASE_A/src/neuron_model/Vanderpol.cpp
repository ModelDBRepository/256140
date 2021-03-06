/***************************************************************************
 *                           Vanderpol.cpp                                 *
 *                           -------------------                           *
 * copyright            : (C) 2013 by Francisco Naveros                    *
 * email                : fnaveros@atc.ugr.es                              *
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "../../include/neuron_model/Vanderpol.h"
#include "../../include/neuron_model/VectorNeuronState.h"

#include <iostream>
#include <cmath>
#include <string>

#include "../../include/openmp/openmp.h"

#include "../../include/spike/EDLUTFileException.h"
#include "../../include/spike/Neuron.h"
#include "../../include/spike/InternalSpike.h"
#include "../../include/spike/PropagatedSpike.h"
#include "../../include/spike/Interconnection.h"

#include "../../include/simulation/Utils.h"

void Vanderpol::LoadNeuronModel(string ConfigFile) throw (EDLUTFileException){
	FILE *fh;
	long Currentline = 0L;
	fh=fopen(ConfigFile.c_str(),"rt");
	if(fh){
		Currentline=1L;
		skip_comments(fh,Currentline);

		this->InitialState = (VectorNeuronState *) new VectorNeuronState(2, true);

		//INTEGRATION METHOD
		this->integrationMethod = LoadIntegrationMethod::loadIntegrationMethod((TimeDrivenNeuronModel *)this, fh, &Currentline,N_NeuronStateVariables, N_DifferentialNeuronState, N_TimeDependentNeuronState);
	}
}

void Vanderpol::SynapsisEffect(int index, VectorNeuronState * State, Interconnection * InputConnection){

}

Vanderpol::Vanderpol(string NeuronTypeID, string NeuronModelID): TimeDrivenNeuronModel(NeuronTypeID, NeuronModelID){
}

Vanderpol::~Vanderpol(void)
{
}

void Vanderpol::LoadNeuronModel() throw (EDLUTFileException){
	this->LoadNeuronModel(this->GetModelID()+".cfg");
}

VectorNeuronState * Vanderpol::InitializeState(){
	return this->GetVectorNeuronState();
}


InternalSpike * Vanderpol::ProcessInputSpike(Interconnection * inter, Neuron * target, double time){
	return 0;
}

	

bool Vanderpol::UpdateState(int index, VectorNeuronState * State, double CurrentTime){
	
	bool * internalSpike=State->getInternalSpike();
	int Size=State->GetSizeState();
	float last_update;
	float elapsed_time;
	float elapsed_time_f;
	int i;
	int CPU_thread_index;

	float * NeuronState;


	for (int i=0; i< Size; i++){

		last_update = State->GetLastUpdateTime(i);
		elapsed_time = CurrentTime - last_update;
		elapsed_time_f=elapsed_time;
		State->AddElapsedTime(i,elapsed_time);

		NeuronState=State->GetStateVariableAt(i);

		this->integrationMethod->NextDifferentialEcuationValue(i, NeuronState, elapsed_time_f);

		State->SetLastUpdateTime(i,CurrentTime);
	}

	return false;

}



ostream & Vanderpol::PrintInfo(ostream & out){
	return out;
}	


void Vanderpol::InitializeStates(int N_neurons, int OpenMPQueueIndex){
	//Initialize neural state variables.
	float initialization[] = {1.085f,-6.0e-3f};
	InitialState->InitializeStates(N_neurons, initialization);

	//Initialize integration method state variables.
	this->integrationMethod->InitializeStates(N_neurons, initialization);
}



void Vanderpol::EvaluateDifferentialEcuation(float * NeuronState, float * AuxNeuronState){
			AuxNeuronState[0]=NeuronState[1];
			AuxNeuronState[1]=1000*(1.0f-NeuronState[0]*NeuronState[0])*NeuronState[1]-NeuronState[0];

}

void Vanderpol::EvaluateTimeDependentEcuation(float * NeuronState, float elapsed_time){
}

int Vanderpol::CheckSynapseTypeNumber(int Type){
	if(Type<N_TimeDependentNeuronState && Type>=0){
		return Type;
	}else{
		cout<<"Neuron model "<<this->GetTypeID()<<", "<<this->GetModelID()<<" does not support input synapses of type "<<Type<<endl;
		return 0;
	}
}

