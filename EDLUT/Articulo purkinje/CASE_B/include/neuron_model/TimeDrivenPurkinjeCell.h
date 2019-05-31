/***************************************************************************
 *                           TimeDrivenPurkinjeCell.h                      *
 *                           -------------------                           *
 * copyright            : (C) 2015 by Richard Carrill, Niceto Luque and    *
						  Francisco Naveros								   *
 * email                : rcarrillo@ugr.es, nluque@ugr.es and			   *
						  fnaveros@atc.ugr.es							   *
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef TIMEDRIVENPURKINJECELL_H_
#define TIMEDRIVENPURKINJECELL_H_

/*!
 * \file TimeDrivenPurkinjeCell.h
 *
 * \author Richard Carrillo
 * \author Niceto Luque
 * \author Francisco Naveros
 * \date May 2015
 *
 * This file declares a class which implement a Purkinje cell model.
 */

#include "./TimeDrivenNeuronModel.h"

#include <string>

using namespace std;

class InputSpike;
class VectorNeuronState;
class Interconnection;



/*!
 * \class TimeDrivenPurkinjeCell.h
 *
 * \brief Time driven neuron model with a membrane potential, five currents and two conductances.
 *
 * \author Richard Carrillo
 * \author Niceto Luque
 * \author Francisco Naveros
 * \date May 2015
 */
class TimeDrivenPurkinjeCell : public TimeDrivenNeuronModel {
	protected:

		static float * alpha_ca_table;
		static float * beta_ca_table;
		static float * alpha_M_table;
		static float * beta_M_table;
		static const float Max_V;
		static const float Min_V;
		static const int TableSize=1024*1024;
		static const float aux;

		/*!
		 * \brief leak current 
		 */
		const float g_L;

		/*!
		 * \brief high-threshold noninactivating calcium current 
		 */
		const float g_Ca;

		/*!
		 * \brief muscarinic receptor suppressed potassium current (or M current) 
		 */
		const float g_M;

		
		/*!
		 * \brief Cylinder length of the soma (cm)
		 */
		const float Cylinder_length_of_the_soma;

				
		/*!
		 * \brief Radius of the soma (cm)
		 */
		const float Radius_of_the_soma;
					
		/*!
		 * \brief Cell area (cm^2)
		 */
		const float Area;
		const float inv_Area;
					
		/*!
		 * \brief Membrane capacitance(uF/cm^2)
		 */
		const float Membrane_capacitance;
		const float inv_Membrane_capacitance;



		/*!
		 * \brief Excitatory reversal potential (mV)
		 */
		float eexc;

		/*!
		 * \brief Inhibitory reversal potential (mV)
		 */
		float einh;

		/*!
		 * \brief Firing threshold (mV)
		 */
		float vthr;

		/*!
		 * \brief Resting potential (mV)
		 */
		float erest;

		/*!
		 * \brief AMPA receptor time constant (ms)
		 */
		float texc;
		float inv_texc;

		/*!
		 * \brief GABA receptor time constant (ms)
		 */
		float tinh;
		float inv_tinh;


		/*!
		 * \brief Refractory period (ms)
		 */
		float tref;
		float tref_0_5;
		float inv_tref_0_5;

		/*!
		 * \brief Peak amplitude (mV)
		 */
		float spkpeak;




		/*!
		 * \brief It loads the neuron model description.
		 *
		 * It loads the neuron type description from the file .cfg.
		 *
		 * \param ConfigFile Name of the neuron description file (*.cfg).
		 *
		 * \throw EDLUTFileException If something wrong has happened in the file load.
		 */
		void LoadNeuronModel(string ConfigFile) throw (EDLUTFileException);


		/*!
		 * \brief It abstracts the effect of an input spike in the cell.
		 *
		 * It abstracts the effect of an input spike in the cell.
		 *
		 * \param index The cell index inside the VectorNeuronState.
		 * \param State Cell current state.
		 * \param InputConnection Input connection from which the input spike has got the cell.
		 */
		void SynapsisEffect(int index, Interconnection * InputConnection);



	public:

		/*!
		 * \brief Number of state variables for each cell.
		*/
		static const int N_NeuronStateVariables=6;

		/*!
		 * \brief Number of state variables witch are calculate with a differential equation for each cell (V, ca and M).
		*/
		static const int N_DifferentialNeuronState=3;

		/*!
		 * \brief Number of state variables witch are calculate with a time dependent equation for each cell (g_exc and g_inh).
		*/
		static const int N_TimeDependentNeuronState=3;


		/*!
		 * \brief Default constructor with parameters.
		 *
		 * It generates a new neuron model object without being initialized.
		 *
		 * \param NeuronTypeID Neuron model identificator.
		 * \param NeuronModelID Neuron model configuration file.
		 */
		TimeDrivenPurkinjeCell(string NeuronTypeID, string NeuronModelID);


		/*!
		 * \brief Class destructor.
		 *
		 * It destroys an object of this class.
		 */
		virtual ~TimeDrivenPurkinjeCell();


		/*!
		 * \brief It loads the neuron model description and tables (if necessary).
		 *
		 * It loads the neuron model description and tables (if necessary).
 		 *
		 * \throw EDLUTFileException If something wrong has happened in the file load.
		 */
		virtual void LoadNeuronModel() throw (EDLUTFileException);


		/*!
		 * \brief It return the Neuron Model VectorNeuronState 
		 *
		 * It return the Neuron Model VectorNeuronState 
		 *
		 */
		virtual VectorNeuronState * InitializeState();


		/*!
		 * \brief It processes a propagated spike (input spike in the cell).
		 *
		 * It processes a propagated spike (input spike in the cell).
		 *
		 * \note This function doesn't generate the next propagated spike. It must be externally done.
		 *
		 * \param inter the interconection which propagate the spike
		 * \param target the neuron which receives the spike
		 * \param time the time of the spike.
		 *
		 * \return A new internal spike if someone is predicted. 0 if none is predicted.
		 */
		virtual InternalSpike * ProcessInputSpike(Interconnection * inter, Neuron * target, double time);


		/*!
		 * \brief Update the neuron state variables.
		 *
		 * It updates the neuron state variables.
		 *
		 * \param index The cell index inside the VectorNeuronState. if index=-1, updating all cell.
		 * \param The current neuron state.
		 * \param CurrentTime Current time.
		 *
		 * \return True if an output spike have been fired. False in other case.
		 */
		virtual bool UpdateState(int index, VectorNeuronState * State, double CurrentTime);



		/*!
		 * \brief It prints the time-driven model info.
		 *
		 * It prints the current time-driven model characteristics.
		 *
		 * \param out The stream where it prints the information.
		 *
		 * \return The stream after the printer.
		 */
		virtual ostream & PrintInfo(ostream & out);


		/*!
		 * \brief It initialice VectorNeuronState.
		 *
		 * It initialice VectorNeuronState.
		 *
		 * \param N_neurons cell number inside the VectorNeuronState.
		 */
		virtual void InitializeStates(int N_neurons, int OpenMPQueueIndex);


		/*!
		 * \brief It evaluates the differential equation in NeuronState and it stores the results in AuxNeuronState.
		 *
		 * It evaluates the differential equation in NeuronState and it stores the results in AuxNeuronState.
		 *
		 * \param NeuronState value of the neuron state variables where differential equations are evaluated.
		 * \param AuxNeuronState results of the differential equations evaluation.
		 */
		virtual void EvaluateDifferentialEcuation(float * NeuronState, float * AuxNeuronState);


		/*!
		 * \brief It evaluates the time depedendent ecuation in NeuronState for elapsed_time and it stores the results in NeuronState.
		 *
		 * It evaluates the time depedendent ecuation in NeuronState for elapsed_time and it stores the results in NeuronState.
		 *
		 * \param NeuronState value of the neuron state variables where time dependent equations are evaluated.
		 * \param elapsed_time integration time step.
		 */
		virtual void EvaluateTimeDependentEcuation(float * NeuronState, float elapsed_time);


		/*!
		 * \brief It Checks if the neuron model has this connection type.
		 *
		 * It Checks if the neuron model has this connection type.
		 *
		 * \param Type input connection type.
		 *
		 * \return A a valid connection type for this neuron model.
		 */
		virtual int CheckSynapseTypeNumber(int Type);




		static float * Generate_alpha_ca_table(){
			float * NewLookUpTable=new float[TableSize];
			for(int i=0; i<TableSize; i++){
				float V = Min_V + ((Max_V-Min_V)*i)/(TableSize-1);
				NewLookUpTable[i]=1.6f/(1+exp(-0.072f*(V-5.0f)));;		
			}
			return NewLookUpTable;
		}

		
		static float * Generate_beta_ca_table(){
			float * NewLookUpTable=new float[TableSize];
			for(int i=0; i<TableSize; i++){
				float V = Min_V + ((Max_V-Min_V)*i)/(TableSize-1);
				NewLookUpTable[i]=(0.02f*(V+8.9f))/(exp((V+8.9f)*0.2f)-1.0f);		
			}
			return NewLookUpTable;
		}

		
		static float * Generate_alpha_M_table(){
			float * NewLookUpTable=new float[TableSize];
			for(int i=0; i<TableSize; i++){
				float V = Min_V + ((Max_V-Min_V)*i)/(TableSize-1);
				NewLookUpTable[i]=0.3f/(1+exp((-V-2.0f)*0.2f));		
			}
			return NewLookUpTable;
		}

		

		static float * Generate_beta_M_table(){
			float * NewLookUpTable=new float[TableSize];
			for(int i=0; i<TableSize; i++){
				float V = Min_V + ((Max_V-Min_V)*i)/(TableSize-1);
				NewLookUpTable[i]=0.001f*exp((-V-60.0f)*0.055555555555555f);		
			}
			return NewLookUpTable;
		}

		static void Get_alpha_and_beta_values(float value, float * alpha_ca, float * beta_ca, float * alpha_M, float * beta_M){
				int position=int((value-Min_V)*aux);
				*alpha_ca=alpha_ca_table[position];
				*beta_ca=beta_ca_table[position];
				*alpha_M=alpha_M_table[position];
				*beta_M=beta_M_table[position];
		} 	

};

#endif /* TIMEDRIVENPURKINJECELL_H_ */
