var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = Network_qse","category":"page"},{"location":"#Network_qse.jl-Documentation","page":"Home","title":"Network_qse.jl Documentation","text":"","category":"section"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Network_qse]","category":"page"},{"location":"#Network_qse.AtomicProperties","page":"Home","title":"Network_qse.AtomicProperties","text":"AtomicProperties\n\nStores all physical properties of an  element with charge number Z and atomic  number Z. Mass in units of MeV with conversion factor: M * [meverg] / c^2 = [erg]\n\n\n\n\n\n","category":"type"},{"location":"#Network_qse.MultiNewtonRaphson-Tuple{Array{T,1} where T,Any,Any,Any,Any}","page":"Home","title":"Network_qse.MultiNewtonRaphson","text":"MultiNewtonRaphson(x::Vector, T, rho, y, ap)\n\n[-7.0010491,-9.1] ∂/∂xᵢ [f₁,...,fₙ] good guess: [-4.394094904641707, -12.915712928215058]\n\n\n\n\n\n","category":"method"},{"location":"#Network_qse.df_nse_condition-NTuple{5,Any}","page":"Home","title":"Network_qse.df_nse_condition","text":"df_nse_condition!(J,μ, T,rho, ap)\n\ncomputes Jacobian ∇f ∈ ℝ²×ℝ² with f ∈ ℝ², μ ∈ ℝ² Contains all 4 partial derivates with respect to the  proton and neutron chemical potentials. \n\n\n\n\n\n","category":"method"},{"location":"#Network_qse.extract_partition_function-Tuple{}","page":"Home","title":"Network_qse.extract_partition_function","text":"extract_partition_function()\n\nchecks indentical charge and atomic numbers in massfrdm and partfrdm. Interpolates ω(T)\n\n\n\n\n\n","category":"method"},{"location":"#Network_qse.initial_guess-Tuple{Any}","page":"Home","title":"Network_qse.initial_guess","text":"initial_guess(ap_ni56)\n\nindex: 438 or 863 or 762 inverse of saha equation with only one species and μₚ = μₙ. returns μ.  We assume that the gas is only composed of Nickel 56 and then solve the Saha Equation for the chemical potentials of Protons and neutrons (which we assume as equal). \n\n\n\n\n\n","category":"method"},{"location":"#Network_qse.logsumexp-Tuple{Any}","page":"Home","title":"Network_qse.logsumexp","text":"logsumexp(arr)\n\nmax + log(sum(exp{arr - max})) instead of log(sum(exp(arr)))\n\n\n\n\n\n","category":"method"},{"location":"#Network_qse.nse_condition-Tuple{Any,Real,Real,Real,Any}","page":"Home","title":"Network_qse.nse_condition","text":"nse_condition!(res, μ, T, ρ, y, ap; precision=400)\n\nreturns zero function, as obtained from condition of  mass conservation and charge neutrality for NSE:\n\nlog (∑ᵢXᵢ) = 0 \nlog(∑ᵢ(Zᵢ/Aᵢ)Xᵢ / y) = 0\n\n\n\n\n\n","category":"method"},{"location":"#Network_qse.prefactor-Tuple{Any}","page":"Home","title":"Network_qse.prefactor","text":"prefactor(pf)\n\nreturns prefactor of Xi, as given in http://cococubed.asu.edu/codepages/nse.shtml The physical quantities to compute theBoltzmann  prefactor  for a given element are temperature,  density (which are varied),  mass and atomic number.\n\n\n\n\n\n","category":"method"},{"location":"#Network_qse.read_mass_frdm-Tuple{}","page":"Home","title":"Network_qse.read_mass_frdm","text":"read_mass_frdm()\n\nread ground state energies\n\n\n\n\n\n","category":"method"},{"location":"#Network_qse.read_part_frdm-Tuple{}","page":"Home","title":"Network_qse.read_part_frdm","text":"read_part_frdm()\n\nread in follwing files:\n\nmass-frdm95.dat:\n\nGround state properties based on the FRDM model Format ––– Each record of the file contains:\n\nZ    : charge number    A    : mass number    El   : element symbol    fl   : flag corresponding to 0 if no experimental data available                                 1 for a mass excess recommended by                                   Audi&Wapstra (1995)                                 2 for a measured mass from                                   Audi&Wapstra (1995)    Mexp : experimental or recommended atomic mass excess in MeV of           Audi&Wapstra (1995)    Mth  : calculated FRDM atomic mass excess in MeV    Emic : calculated FRDM microscopic energy in MeV    beta2: calculated quadrupole deformation of the nuclear ground-state    beta3: calculated octupole deformation of the nuclear ground-state    beta4: calculated hexadecapole deformation of the nuclear ground-state    beta6: calculated hexacontatetrapole deformation of the nuclear           ground-state\n\nThe corresponding FORTRAN format is (2i4,1x,a2,1x,i1,3f10.3,4f8.3)\n\n\n\npart_frdm.asc:\n\nEach of the two files starts with 4 header lines briefly summarizing the contents of the given columns. This is followed by entries sorted by charge and mass number of the isotope. Each table ends with the line \"END OF TABLE\". Each entry consists of 5 lines:\n\nIsotope (in standard notation);\nCharge number of isotope, mass number of isotope, ground state spin of the isotope;\n\n3-5. Partition functions normalized to the g.s. spin;    Third line: Partition functions for the  temperatures (in 10^9 K):    0.01, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7;        Fourth line:  Partition functions for the temperatures (in 10^9 K):    0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5;        Fifth line: Partition functions for the temperatures (in 10^9 K):    4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 . Information for the next isotope starts after the last partition function line.\n\nPIA: ADDED partition functions and A,Z,s for n,p,d,T,3He\n\n\n\n\n\n","category":"method"},{"location":"#Network_qse.read_species-Tuple{}","page":"Home","title":"Network_qse.read_species","text":"read_species()\n\nreads species from table\n\n\n\n\n\n","category":"method"},{"location":"#Notes","page":"Home","title":"Notes","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This layout should be customized in docs/src/index.md","category":"page"}]
}
