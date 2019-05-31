/***************************************************************************
 *                           LIFTimeDrivenModel_1_2.cpp                    *
 *                           -------------------                           *
 * copyright            : (C) 2013 by Jesus Garrido and Francisco Naveros  *
 * email                : jgarrido@atc.ugr.es, fnaveros@atc.ugr.es         *
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "../../include/neuron_model/TimeDrivenPurkinjeCell.h"
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

#include "../../include/openmp/openmp.h"


void TimeDrivenPurkinjeCell::LoadNeuronModel(string ConfigFile) throw (EDLUTFileException){
	FILE *fh;
	long Currentline = 0L;
	fh=fopen(ConfigFile.c_str(),"rt");
	if(fh){
		//Currentline=1L;
		//skip_comments(fh,Currentline);
		//if(fscanf(fh,"%f",&this->eexc)==1){
		//	skip_comments(fh,Currentline);

		//	if (fscanf(fh,"%f",&this->einh)==1){
		//		skip_comments(fh,Currentline);

		//		if(fscanf(fh,"%f",&this->erest)==1){
		//			skip_comments(fh,Currentline);

		//			if(fscanf(fh,"%f",&this->vthr)==1){
		//				skip_comments(fh,Currentline);

		//				if(fscanf(fh,"%f",&this->cm)==1){
		//					inv_cm=1.0f/cm;
		//					skip_comments(fh,Currentline);

		//					if(fscanf(fh,"%f",&this->texc)==1){
		//						inv_texc=1.0f/texc;
		//						skip_comments(fh,Currentline);

		//						if(fscanf(fh,"%f",&this->tinh)==1){
		//							inv_tinh=1.0f/tinh;
		//							skip_comments(fh,Currentline);

		//							if(fscanf(fh,"%f",&this->tref)==1){
		//								skip_comments(fh,Currentline);

		//								if(fscanf(fh,"%f",&this->grest)==1){
		//									skip_comments(fh,Currentline);

											this->InitialState = (VectorNeuronState *) new VectorNeuronState(N_NeuronStateVariables, true);

		//								} else {
		//									throw EDLUTFileException(13,60,3,1,Currentline);
		//								}
		//							} else {
		//								throw EDLUTFileException(13,61,3,1,Currentline);
		//							}
		//						} else {
		//							throw EDLUTFileException(13,62,3,1,Currentline);
		//						}
		//					} else {
		//						throw EDLUTFileException(13,63,3,1,Currentline);
		//					}
		//				} else {
		//					throw EDLUTFileException(13,64,3,1,Currentline);
		//				}
		//			} else {
		//				throw EDLUTFileException(13,65,3,1,Currentline);
		//			}
		//		} else {
		//			throw EDLUTFileException(13,66,3,1,Currentline);
		//		}
		//	} else {
		//		throw EDLUTFileException(13,67,3,1,Currentline);
		//	}
		//} else {
		//	throw EDLUTFileException(13,68,3,1,Currentline);
		//}

		//INTEGRATION METHOD
		this->integrationMethod = LoadIntegrationMethod::loadIntegrationMethod((TimeDrivenNeuronModel *)this, fh, &Currentline, N_NeuronStateVariables, N_DifferentialNeuronState, N_TimeDependentNeuronState);

	}
}

void TimeDrivenPurkinjeCell::SynapsisEffect(int index, Interconnection * InputConnection){
	this->GetVectorNeuronState()->IncrementStateVariableAtCPU(index,N_DifferentialNeuronState+InputConnection->GetType(),1e-6f*InputConnection->GetWeight());
}

TimeDrivenPurkinjeCell::TimeDrivenPurkinjeCell(string NeuronTypeID, string NeuronModelID): TimeDrivenNeuronModel(NeuronTypeID, NeuronModelID), g_L(0.02f),
		g_Ca(0.001f), g_M(0.75f), Cylinder_length_of_the_soma(0.0015f), Radius_of_the_soma(0.0008f), Area(3.141592f*0.0015f*2.0f*0.0008f),
		inv_Area(1.0f/(3.141592f*0.0015f*2.0f*0.0008f)), Membrane_capacitance(1.0f), inv_Membrane_capacitance(1.0f/1.0f)
		{
	eexc=0.0f;
	einh=-80.0f ;
	vthr= -35.0f;
	erest=-65.0f;
	texc=1.0f;
	inv_texc=1.0f/texc;
	tinh=2;
	inv_tinh=1.0f/tinh;
	tref=1.35f;
	tref_0_5=tref*0.5f;
	inv_tref_0_5=1.0f/tref_0_5;
	spkpeak=31.0f;
}

TimeDrivenPurkinjeCell::~TimeDrivenPurkinjeCell(void)
{
}

void TimeDrivenPurkinjeCell::LoadNeuronModel() throw (EDLUTFileException){
	this->LoadNeuronModel(this->GetModelID()+".cfg");
}

VectorNeuronState * TimeDrivenPurkinjeCell::InitializeState(){
	return this->GetVectorNeuronState();
}


InternalSpike * TimeDrivenPurkinjeCell::ProcessInputSpike(Interconnection * inter, Neuron * target, double time){
	// Add the effect of the input spike
	this->SynapsisEffect(target->GetIndex_VectorNeuronState(),inter);

	return 0;
}






bool TimeDrivenPurkinjeCell::UpdateState(int index, VectorNeuronState * State, double CurrentTime){
	//float * NeuronState;
	//NeuronState[0] --> V
	//NeuronState[1] --> c
	//NeuronState[2] --> M
	//NeuronState[3] --> gexc 
	//NeuronState[4] --> ginh 


	bool * internalSpike=State->getInternalSpike();

	double last_update = State->GetLastUpdateTime(0);
	double elapsed_time = CurrentTime - last_update;
	float elapsed_time_f=elapsed_time;


	for(int j=0; j<NumberOfOpenMPTasks-1; j++){
		#ifdef _OPENMP 
			#if	_OPENMP >= OPENMPVERSION30
				#pragma omp task firstprivate (j) shared(internalSpike, State, CurrentTime) 
			#endif
		#endif
		{
			for (int i=LimitOfOpenMPTasks[j]; i< LimitOfOpenMPTasks[j+1]; i++){
				State->AddElapsedTime(i,elapsed_time);
				double last_spike = State->GetLastSpikeTime(i);
				last_spike*=1000;//ms
		 
				float * NeuronState=State->GetStateVariableAt(i);
					
				bool spike = false;
NeuronState[5]=last_spike;
				this->integrationMethod->NextDifferentialEcuationValue(i, NeuronState, elapsed_time_f*1000);



				if (last_spike > this->tref && NeuronState[0] > this->vthr) {
					State->NewFiredSpike(i);
					last_spike=0;
					spike = true;
				}


				if(last_spike < tref){
					if(last_spike <= tref_0_5){
						NeuronState[0]=vthr+(spkpeak-vthr)*(last_spike*inv_tref_0_5);
					}else{
						NeuronState[0]=spkpeak-(spkpeak-erest)*((last_spike-tref_0_5)*inv_tref_0_5);
					}
				}

				internalSpike[i]=spike;
				State->SetLastUpdateTime(i,CurrentTime);
			}
		}
	}
	
	for (int i=LimitOfOpenMPTasks[NumberOfOpenMPTasks-1]; i< LimitOfOpenMPTasks[NumberOfOpenMPTasks]; i++){
		State->AddElapsedTime(i,elapsed_time);
		double last_spike = State->GetLastSpikeTime(i);
		last_spike*=1000;//ms
 
		float * NeuronState=State->GetStateVariableAt(i);
			
		bool spike = false;
NeuronState[5]=last_spike;
		this->integrationMethod->NextDifferentialEcuationValue(i, NeuronState, elapsed_time_f*1000);



		if (last_spike > this->tref && NeuronState[0] > this->vthr) {
			State->NewFiredSpike(i);
			last_spike=0;
			spike = true;
		}

		if(last_spike < tref){
			if(last_spike <= tref_0_5){
				NeuronState[0]=vthr+(spkpeak-vthr)*(last_spike*inv_tref_0_5);
			}else{
				NeuronState[0]=spkpeak-(spkpeak-erest)*((last_spike-tref_0_5)*inv_tref_0_5);
			}
		}


		internalSpike[i]=spike;

		State->SetLastUpdateTime(i,CurrentTime);
	}

	#ifdef _OPENMP 
		#if	_OPENMP >= OPENMPVERSION30
			#pragma omp taskwait
		#endif
	#endif

	return false;
}



ostream & TimeDrivenPurkinjeCell::PrintInfo(ostream & out){
	//out << "- Leaky Time-Driven Model: " << this->GetModelID() << endl;

	//out << "\tExc. Reversal Potential: " << this->eexc << "V\tInh. Reversal Potential: " << this->einh << "V\tResting potential: " << this->erest << "V" << endl;

	//out << "\tFiring threshold: " << this->vthr << "V\tMembrane capacitance: " << this->cm << "nS\tExcitatory Time Constant: " << this->texc << "s" << endl;

	//out << "\tInhibitory time constant: " << this->tinh << "s\tRefractory Period: " << this->tref << "s\tResting Conductance: " << this->grest << "nS" << endl;

	return out;
}	



void TimeDrivenPurkinjeCell::InitializeStates(int N_neurons, int OpenMPQueueIndex){
	//Initialize neural state variables.
	float alpha_ca, beta_ca, alpha_M, beta_M;
	Get_alpha_and_beta_values(erest, &alpha_ca, &beta_ca, &alpha_M, &beta_M);
	//c_inf
	float c_inf=alpha_ca/(alpha_ca+beta_ca);

	//M_inf
	float M_inf=alpha_M/(alpha_M+beta_M);

	float initialization[] = {erest,c_inf,M_inf,0.0f,0.0f};
	InitialState->InitializeStates(N_neurons, initialization);

	//Initialize integration method state variables.
	this->integrationMethod->InitializeStates(N_neurons, initialization);


	//Calculate number of OpenMP task and size of each one.
	CalculateTaskSizes(N_neurons, 50);
}



void TimeDrivenPurkinjeCell::EvaluateDifferentialEcuation(float * NeuronState, float * AuxNeuronState){
	float V=NeuronState[0];
	float ca=NeuronState[1];
	float M=NeuronState[2];
	float g_exc=NeuronState[3];
	float g_inh=NeuronState[4];
	float last_spike=NeuronState[5];

	//V
	if(last_spike >= tref){
		AuxNeuronState[0]=(-g_L*(V+70.0f)-g_Ca*ca*ca*(V-125.0f)-g_M*M*(V+95.0f) + (g_exc * (this->eexc - V) + g_inh * (this->einh - V))*inv_Area )*inv_Membrane_capacitance;
	}else if(last_spike <= tref_0_5){
		AuxNeuronState[0]=(spkpeak-vthr)*inv_tref_0_5;
	}else{
		AuxNeuronState[0]=(erest-spkpeak)*inv_tref_0_5;
	}

	float alpha_ca, beta_ca, alpha_M, beta_M;
	Get_alpha_and_beta_values(V, &alpha_ca, &beta_ca, &alpha_M, &beta_M);

	//ca
	float inv_tau_ca=(alpha_ca+beta_ca);
	AuxNeuronState[1]=alpha_ca - ca*inv_tau_ca;
	
	//M	
	float inv_tau_M=(alpha_M+beta_M);
	AuxNeuronState[2]=alpha_M - M*inv_tau_M;


}

void TimeDrivenPurkinjeCell::EvaluateTimeDependentEcuation(float * NeuronState, float elapsed_time){
	float limit=1e-20;
	
	if(NeuronState[N_DifferentialNeuronState]<limit){
		NeuronState[N_DifferentialNeuronState]=0.0f;
	}else{
		NeuronState[N_DifferentialNeuronState]*= ExponentialTable::GetResult(-(elapsed_time*this->inv_texc));
	}
	if(NeuronState[N_DifferentialNeuronState+1]<limit){
		NeuronState[N_DifferentialNeuronState+1]=0.0f;
	}else{
		NeuronState[N_DifferentialNeuronState+1]*= ExponentialTable::GetResult(-(elapsed_time*this->inv_tinh));
	}	
}


int TimeDrivenPurkinjeCell::CheckSynapseTypeNumber(int Type){
	if(Type<N_TimeDependentNeuronState && Type>=0){
		return Type;
	}else{
		cout<<"Neuron model "<<this->GetTypeID()<<", "<<this->GetModelID()<<" does not support input synapses of type "<<Type<<endl;
		return 0;
	}
}



const float TimeDrivenPurkinjeCell::Max_V=35.0f;
const float TimeDrivenPurkinjeCell::Min_V=-100.0f;

const float TimeDrivenPurkinjeCell::aux=(TimeDrivenPurkinjeCell::TableSize-1)/( TimeDrivenPurkinjeCell::Max_V- TimeDrivenPurkinjeCell::Min_V);

float * TimeDrivenPurkinjeCell::alpha_ca_table=Generate_alpha_ca_table();
float * TimeDrivenPurkinjeCell::beta_ca_table=Generate_beta_ca_table();
float * TimeDrivenPurkinjeCell::alpha_M_table=Generate_alpha_M_table();
float * TimeDrivenPurkinjeCell::beta_M_table=Generate_beta_M_table();


