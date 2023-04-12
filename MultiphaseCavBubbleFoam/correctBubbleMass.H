 if (correctBubbleMass)
{
    forAllConstIters(mixture.phases(), phase)
    {
        if (phase().name() == "air")
        {
            alpha2rho = phase()*phase().thermo().rho();
            bubblem = Foam::gSum(mesh.V()*alpha2rho.field());
            bubbleV = Foam::gSum(mesh.V()*phase().field());
            for (phaseModel& PHASE : mixture.phases())
            {
                if (PHASE.name() == "air")
                {
                    PHASE.thermo().correctRho(phase().thermo().rho()*(bubblem0/bubblem - 1), rhoMin, rhoMax);
                }
            }
            if (bubblePhase == 0)
            {
                if (bubbleV < bubbleV0)
                {
                    bubblePhase = 1;
                    for (phaseModel& PHASE : mixture.phases())
                    {
                        if (PHASE.name() == "air")
                        {
                            PHASE.thermo().correctRho(phase().thermo().rho()*(pFact - 1), rhoMin, rhoMax);
                        }
                    }
                    alpha2rho = phase()*phase().thermo().rho();
                    bubblem0 = Foam::gSum(mesh.V()*alpha2rho.field());
                }
                bubbleV0 = bubbleV;
            } else if (bubblePhase == 1)
            {
                if (bubbleV > bubbleV0)
                {
                    bubblePhase = 2;
                }
                bubbleV0 = bubbleV;
            }
            p_g = max((p0_air + B_air) * pow(phase().thermo().rho()/rho0_air , gamma_air ) - B_air,pMin);
            p_rgh = (1 - phase()) * p_rgh + phase() * p_g;
            p = max(p_rgh + mixture.rho()*gh, pMin);
            for (phaseModel& phase : mixture.phases())
            {
                if (phase.name() == "water")
                {
                    mixture.correctComponent(phase);
                }
            }
    
            forAllConstIters(mixture.phases(), phase)
            {
                if (phase().name() == "water")
                {
                    p_temp = p*phase();
                }
            }
    
            p = max(p,pMinGas);
    
            forAllConstIters(mixture.phases(), phase)
            {
                if (phase().name() == "water")
                {
                    p_temp = p_temp + p*(1-phase());
                }
            }
    
            for (phaseModel& phase : mixture.phases())
            {
                if (phase.name() == "oil")
                {
                    mixture.correctComponent(phase);
                }
                else if (phase.name() == "air")
                {
                    mixture.correctComponent(phase);
                }
            }

            mixture.correctFinal();

            p = p_temp;
        }
    }
}
