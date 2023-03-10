#include "HepMC3/GenVertex.h"
#include "HepMC3/GenVertex_fwd.h"
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"

#include "inputParameters.h"
#include "upcXevent.h"
#include "hepmc3writer.h"

using HepMC3::FourVector;
using HepMC3::Print;
using namespace starlightConstants;
using namespace std;

hepMC3Writer::hepMC3Writer()
{}
int hepMC3Writer::initWriter(const inputParameters &param)
{
    std::string hepmc3_filename = param.baseFileName() + ".hepmc";

    //setting imporant collision parameters
    initBeamHepMC3(param);
    _hepmc3_output = new HepMC3::WriterAscii(hepmc3_filename);
    return 0;
}

particleTypeEnum getVMPid(particleTypeEnum pid, decayTypeEnum pdecay){
    if(pdecay == LEPTONPAIR){
        return UNKNOWN;//lepton pairs have no intervening Vector meson.
    }
    

    switch(pid){
        case ZOVERZ03:
        case RHO_ee:
        case RHO_mumu:
        case RHOZEUS:
        case RHO:
            return RHO;
        break;
        case TAUONDECAY:
        case TAUON:
            return TAUON;
        break;
        case OMEGA_pipipi:
        case OMEGA:
            return OMEGA;
        break;
        case PHI_ee:
        case PHI:
            return PHI;
        break;
        case JPSI2S_mumu:
        case JPSI2S_ee:
        case JPSI2S:
            return JPSI2S;
        break;
        case JPSI_mumu:
        case JPSI_ee:
        case JPSI:
            return JPSI;
        break;
        case UPSILON_mumu:
        case UPSILON_ee:
        case UPSILON:
            return UPSILON;
        break;
        case UPSILON2S_mumu:
        case UPSILON2S_ee:
        case UPSILON2S:
            return UPSILON2S;
        break;
        default:
            assert(true);

    }
    return pid;
}

static int getNucleusPdgCode(int Z, int A){
    if(A == 1){
        if(Z ==1){
            return 2212;//proton //https://www.star.bnl.gov/public/comp/simu/newsite/gstar/kumacs/NewParticle.html  
        }
        else if (Z==0){
            return 2112;//neutron //https://www.star.bnl.gov/public/comp/simu/newsite/gstar/kumacs/NewParticle.html  
        }
        else{ printWarn << "Error, unkown beam"<< endl;}
    }
    else{
        
        if (Z< 0 || A<1) printWarn << "Error, Unkown beam"<<std::endl;
        
        //https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
        if(A < 1000 && A>Z){
            return 1000000000 + Z* 10000 + A* 10;//Form 10LZZZAAAI //I=0 L=0 for ground state nucleus with no strange quarks.
        }else printWarn << "Impossible Beam"<<std::endl;
    }
    return -1;
}

int hepMC3Writer::initBeamHepMC3(const inputParameters &param)
{
    double protonMass = param.protonMass();
    double beam1E = protonMass * param.beam1A() * param.beam1LorentzGamma();
    double beam2E = protonMass * param.beam2A() * param.beam2LorentzGamma();
    double beam1pz = protonMass* param.beam1A() * sqrt(param.beam1LorentzGamma()*param.beam1LorentzGamma() - 1);
    double beam2pz = - protonMass* param.beam2A() * sqrt(param.beam2LorentzGamma()*param.beam2LorentzGamma() - 1);

    hepmc3_beam1_four_vector.set(0,0,beam1pz,beam1E);
    hepmc3_beam2_four_vector.set(0,0,beam2pz,beam2E);
    
    beam1_pdg_id = getNucleusPdgCode(param.beam1Z(), param.beam1A());
    beam2_pdg_id = getNucleusPdgCode(param.beam2Z(), param.beam2A());
    VM_pdg_id = getVMPid(param.prodParticleType(), param.prodParticleDecayType());
    PID = param.prodParticleType();
    _decay = param.prodParticleDecayType();

    return 0;
}

int hepMC3Writer::writeEvent(const upcXEvent &event, int eventnumber){

    HepMC3::GenEvent hepmc3_evt(HepMC3::Units::GEV, HepMC3::Units::MM);
    hepmc3_evt.set_event_number(eventnumber);

    HepMC3::GenParticlePtr hepmc3_beam1_in = std::make_shared<HepMC3::GenParticle>(hepmc3_beam1_four_vector,beam1_pdg_id,4);
    HepMC3::GenParticlePtr hepmc3_beam2_in = std::make_shared<HepMC3::GenParticle>(hepmc3_beam2_four_vector,beam2_pdg_id,4);

    
    
    if (_decay == NARROWVMDEFAULT || _decay == WIDEVMDEFAULT){
                
        lorentzVector gamma = (*event.getGamma())[0];
        //we can use the status code 13 to represent the virtual photons.//insight drawn from estarlight.
        HepMC3::GenParticlePtr gamma_particle = std::make_shared<HepMC3::GenParticle>(FourVector(gamma.GetPx(),
                                                                                            gamma. GetPy(),
                                                                                            gamma.GetPz(),
                                                                                            gamma.GetE()),PHOTON,13);

        
        lorentzVector beam_1 = event.getBeam1();
        lorentzVector beam_2 = event.getBeam2();
        HepMC3::GenParticlePtr hepmc3_beam1_out = std::make_shared<HepMC3::GenParticle>(FourVector(beam_1.GetPx(),
                                                                                                    beam_1.GetPy(),
                                                                                                    beam_1.GetPz(),
                                                                                                    beam_1.GetE()),beam1_pdg_id,1);
        HepMC3::GenParticlePtr hepmc3_beam2_out = std::make_shared<HepMC3::GenParticle>(FourVector(beam_2.GetPx(),
                                                                                                    beam_2.GetPy(),
                                                                                                    beam_2.GetPz(),
                                                                                                    beam_2.GetE()),beam2_pdg_id,1);


        lorentzVector vmeson = (*event.getVectorMeson())[0];
        HepMC3::GenParticlePtr hepmc3_vector_meson = std::make_shared<HepMC3::GenParticle>(FourVector(vmeson.GetPx(),vmeson.GetPy(),vmeson.GetPz(),vmeson.GetE()),VM_pdg_id,2);

        const std::vector<starlightParticle> * particle_vector = event.getParticles();

        HepMC3::GenVertexPtr hepmc3_root_vertex = std::make_shared<HepMC3::GenVertex>(FourVector(0,0,0,0));
        HepMC3::GenVertexPtr hepmc3_gamEmit_vertex = std::make_shared<HepMC3::GenVertex>(FourVector(0,0,0,0));
        HepMC3::GenVertexPtr hepmc3_VMProd_vertex = std::make_shared<HepMC3::GenVertex>(FourVector(0,0,0,0));
        HepMC3::GenVertexPtr hepmc3_PartProd_vertex = std::make_shared<HepMC3::GenVertex>(FourVector(0,0,0,0));

        hepmc3_evt.add_particle(hepmc3_beam1_in);
        hepmc3_evt.add_particle(hepmc3_beam2_in);
        if(event.targetBeamNo() == 1){
            hepmc3_evt.add_particle(gamma_particle);
            hepmc3_evt.add_particle(hepmc3_beam2_out);
            hepmc3_evt.add_particle(hepmc3_beam1_out);
            hepmc3_evt.add_particle(hepmc3_vector_meson);
            hepmc3_evt.add_vertex(hepmc3_root_vertex);
            hepmc3_evt.add_vertex(hepmc3_gamEmit_vertex);
            hepmc3_evt.add_vertex(hepmc3_VMProd_vertex);
            hepmc3_evt.add_vertex(hepmc3_PartProd_vertex);
            
            hepmc3_root_vertex->add_particle_out(hepmc3_beam1_in);
            hepmc3_root_vertex->add_particle_out(hepmc3_beam2_in);

            hepmc3_gamEmit_vertex->add_particle_in(hepmc3_beam2_in);
            hepmc3_gamEmit_vertex->add_particle_out(gamma_particle);
            hepmc3_gamEmit_vertex->add_particle_out(hepmc3_beam2_out);

            hepmc3_VMProd_vertex->add_particle_in(hepmc3_beam1_in);
            hepmc3_VMProd_vertex->add_particle_in(gamma_particle);            
            hepmc3_VMProd_vertex->add_particle_out(hepmc3_beam1_out);
            hepmc3_VMProd_vertex->add_particle_out(hepmc3_vector_meson);

            hepmc3_PartProd_vertex->add_particle_in(hepmc3_vector_meson);


        }else if(event.targetBeamNo() == 2){
            hepmc3_evt.add_particle(hepmc3_beam1_out);
            hepmc3_evt.add_particle(gamma_particle);
            hepmc3_evt.add_particle(hepmc3_vector_meson);
            hepmc3_evt.add_particle(hepmc3_beam2_out);
            hepmc3_evt.add_vertex(hepmc3_root_vertex);
            hepmc3_evt.add_vertex(hepmc3_gamEmit_vertex);
            hepmc3_evt.add_vertex(hepmc3_VMProd_vertex);
            hepmc3_evt.add_vertex(hepmc3_PartProd_vertex);

            hepmc3_root_vertex->add_particle_out(hepmc3_beam1_in);
            hepmc3_root_vertex->add_particle_out(hepmc3_beam2_in);

            hepmc3_gamEmit_vertex->add_particle_in(hepmc3_beam1_in);
            hepmc3_gamEmit_vertex->add_particle_out(hepmc3_beam1_out);
            hepmc3_gamEmit_vertex->add_particle_out(gamma_particle);

            hepmc3_VMProd_vertex->add_particle_in(gamma_particle);
            hepmc3_VMProd_vertex->add_particle_in(hepmc3_beam2_in);
            hepmc3_VMProd_vertex->add_particle_out(hepmc3_vector_meson);                      
            hepmc3_VMProd_vertex->add_particle_out(hepmc3_beam2_out);

            hepmc3_PartProd_vertex->add_particle_in(hepmc3_vector_meson);

        }
        else assert(false);

        for ( std::vector<starlightParticle>::const_iterator particle_iter = (*particle_vector).begin(); particle_iter != (*particle_vector).end(); 	++particle_iter)
        {
            int hepmc3_pid = (*particle_iter).getPdgCode();
            /** pass to HepMC3 FourVector **/
            FourVector hepmc3_four_vector = FourVector( (*particle_iter).GetPx(),
                                (*particle_iter).GetPy(),
                                (*particle_iter).GetPz(),
                                (*particle_iter).GetE());
            
            HepMC3::GenParticlePtr hepmc3_particle = std::make_shared<HepMC3::GenParticle>( hepmc3_four_vector, hepmc3_pid, 1 );
            hepmc3_evt.add_particle(hepmc3_particle);
            hepmc3_PartProd_vertex->add_particle_out(hepmc3_particle);
            
        }
        
    }else if(_decay == LEPTONPAIR || _decay == SINGLEMESON ){
        //1. create all particles

        bool leptonpair = false;//Is this is a two photon to lepton pair interaction?
        if(_decay == LEPTONPAIR) leptonpair = true;

        lorentzVector beam_1 = event.getBeam1();
        HepMC3::GenParticlePtr hepmc3_beam1_out = std::make_shared<HepMC3::GenParticle>(FourVector(beam_1.GetPx(),
                                                                                                    beam_1.GetPy(),
                                                                                                    beam_1.GetPz(),
                                                                                                    beam_1.GetE()),beam1_pdg_id,1);
        

        lorentzVector gamma1 = (*event.getGamma())[0];
        lorentzVector gamma2 = (*event.getGamma())[1];
        //we can use the status code 13 to represent the virtual photons.//insight drawn from estarlight.
        HepMC3::GenParticlePtr gamma_particle1 = std::make_shared<HepMC3::GenParticle>(FourVector(gamma1.GetPx(),
                                                                                                gamma1. GetPy(),
                                                                                                gamma1.GetPz(),
                                                                                                gamma1.GetE()),PHOTON,13);
        HepMC3::GenParticlePtr gamma_particle2 = std::make_shared<HepMC3::GenParticle>(FourVector(gamma2.GetPx(),
                                                                                                gamma2. GetPy(),
                                                                                                gamma2.GetPz(),
                                                                                                gamma2.GetE()),PHOTON,13);
        lorentzVector beam_2 = event.getBeam2();
        HepMC3::GenParticlePtr hepmc3_beam2_out = std::make_shared<HepMC3::GenParticle>(FourVector(beam_2.GetPx(),
                                                                                                    beam_2.GetPy(),
                                                                                                    beam_2.GetPz(),
                                                                                                    beam_2.GetE()),beam2_pdg_id,1);
        
        HepMC3::GenParticlePtr hepmc3_vector_meson;
        if(!leptonpair){
            lorentzVector vmeson = (*event.getVectorMeson())[0];
            hepmc3_vector_meson = std::make_shared<HepMC3::GenParticle>(FourVector(vmeson.GetPx(),vmeson.GetPy(),vmeson.GetPz(),vmeson.GetE()),VM_pdg_id,2);
        }

        const std::vector<starlightParticle> * particle_vector = event.getParticles();
        //decay product particles created later

        //2. create vertex
        HepMC3::GenVertexPtr hepmc3_root_vertex = std::make_shared<HepMC3::GenVertex>(FourVector(0,0,0,0));
        HepMC3::GenVertexPtr hepmc3_beam1gamEmit_vertex = std::make_shared<HepMC3::GenVertex>(FourVector(0,0,0,0));
        HepMC3::GenVertexPtr hepmc3_beam2gamEmit_vertex = std::make_shared<HepMC3::GenVertex>(FourVector(0,0,0,0));

        HepMC3::GenVertexPtr hepmc3_VMProd_vertex;
        if(!leptonpair) hepmc3_VMProd_vertex = std::make_shared<HepMC3::GenVertex>(FourVector(0,0,0,0));

        HepMC3::GenVertexPtr hepmc3_PartProd_vertex = std::make_shared<HepMC3::GenVertex>(FourVector(0,0,0,0));

        //3. add particles to event

        hepmc3_evt.add_particle(hepmc3_beam1_in);
        hepmc3_evt.add_particle(hepmc3_beam2_in);
        hepmc3_evt.add_particle(hepmc3_beam1_out);
        hepmc3_evt.add_particle(gamma_particle1);
        hepmc3_evt.add_particle(gamma_particle2);
        hepmc3_evt.add_particle(hepmc3_beam2_out);
        if(!leptonpair) hepmc3_evt.add_particle(hepmc3_vector_meson);
        //decay product particles added later

        //4. add vertex to event
        hepmc3_evt.add_vertex(hepmc3_root_vertex);
        hepmc3_evt.add_vertex(hepmc3_beam1gamEmit_vertex);
        hepmc3_evt.add_vertex(hepmc3_beam2gamEmit_vertex);
        if(!leptonpair) hepmc3_evt.add_vertex(hepmc3_VMProd_vertex);
        hepmc3_evt.add_vertex(hepmc3_PartProd_vertex);


        //5. add particles to vertex
        hepmc3_root_vertex->add_particle_out(hepmc3_beam1_in);
        hepmc3_root_vertex->add_particle_out(hepmc3_beam2_in);

        hepmc3_beam1gamEmit_vertex->add_particle_in(hepmc3_beam1_in);
        hepmc3_beam1gamEmit_vertex->add_particle_out(hepmc3_beam1_out);
        hepmc3_beam1gamEmit_vertex->add_particle_out(gamma_particle1);

        hepmc3_beam2gamEmit_vertex->add_particle_in(hepmc3_beam2_in);
        hepmc3_beam2gamEmit_vertex->add_particle_out(gamma_particle2);
        hepmc3_beam2gamEmit_vertex->add_particle_out(hepmc3_beam2_out);

        if(!leptonpair){
            hepmc3_VMProd_vertex->add_particle_in(gamma_particle1);
            hepmc3_VMProd_vertex->add_particle_in(gamma_particle2);
            hepmc3_VMProd_vertex->add_particle_out(hepmc3_vector_meson);

            hepmc3_PartProd_vertex->add_particle_in(hepmc3_vector_meson);
        }else{
            hepmc3_PartProd_vertex->add_particle_in(gamma_particle1);
            hepmc3_PartProd_vertex->add_particle_in(gamma_particle2);
        }//product particles handled next.

        //6. dealing with the output produced particles
        for ( std::vector<starlightParticle>::const_iterator particle_iter = (*particle_vector).begin(); particle_iter != (*particle_vector).end(); 	++particle_iter)
        {
            int hepmc3_pid = (*particle_iter).getPdgCode();
            /** pass to HepMC3 FourVector **/
            FourVector hepmc3_four_vector = FourVector( (*particle_iter).GetPx(),
                                (*particle_iter).GetPy(),
                                (*particle_iter).GetPz(),
                                (*particle_iter).GetE());
            
            HepMC3::GenParticlePtr hepmc3_particle = std::make_shared<HepMC3::GenParticle>( hepmc3_four_vector, hepmc3_pid, 1 );//1.particle created
            hepmc3_evt.add_particle(hepmc3_particle);//3.particle added to event
            hepmc3_PartProd_vertex->add_particle_out(hepmc3_particle);//5.particle added to vertex
        }       


    }

    _hepmc3_output->write_event(hepmc3_evt);
    hepmc3_evt.clear();

    return eventnumber;
}

