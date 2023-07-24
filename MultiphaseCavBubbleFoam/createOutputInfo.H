volScalarField alpha2p = p_rgh;
volScalarField alpha1rho = rho;
volScalarField alpha2rho = rho;
volScalarField alpha3rho = rho;
scalar bubblem0 = 0;
scalar bubblem = 0;
scalar bubbleV0 = 0;
scalar bubbleV = 0;
int bubblePhase = 0;

std::ofstream infovstime;
infovstime.open ("info.csv", std::ios::out|std::ios::trunc);
forAllConstIters(mixture.phases(), phase)
{
    if (phase().name() == "oil")
    {
        alpha3rho = phase()*phase().thermo().rho();
    }
    if (phase().name() == "water")
    {
        alpha1rho = phase()*phase().thermo().rho();
    }
}
forAllConstIters(mixture.phases(), phase)
{
    if (phase().name() == "air")
    {
        alpha2p = phase()*p_rgh;
        alpha2rho = phase()*phase().thermo().rho();
        bubbleV0 = Foam::gSum(mesh.V()*phase().field());
        bubblem0 = Foam::gSum(mesh.V()*alpha2rho.field());
        //               1         2            3            4         5        6
        infovstime << "%Time,bubbleVolume,bubblepressure,bubblemass,oilmass,watermass" << nl;
        infovstime << runTime.timeName()<<","
                   << bubbleV0 << ","
                   << Foam::gSum(mesh.V()*alpha2p.field())/bubbleV0 << ","
                   << bubblem0 << ","
                   << Foam::gSum(mesh.V()*alpha3rho.field()) << ","
                   << Foam::gSum(mesh.V()*alpha1rho.field()) << nl;
        infovstime.flush();
    }
}
