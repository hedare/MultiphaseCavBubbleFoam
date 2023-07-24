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
        if (correctBubbleMass == false)
        {
            alpha2rho = phase()*phase().thermo().rho();
            bubblem = Foam::gSum(mesh.V()*alpha2rho.field());
            bubbleV = Foam::gSum(mesh.V()*phase().field());
        }
        infovstime << runTime.timeName()<<","
                   << bubbleV << ","
                   << Foam::gSum(mesh.V()*alpha2p.field())/bubbleV << ","
                   << bubblem << ","
                   << Foam::gSum(mesh.V()*alpha3rho.field()) << ","
                   << Foam::gSum(mesh.V()*alpha1rho.field()) << nl;
        infovstime.flush();
    }
}
