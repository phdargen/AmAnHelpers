class Particle:
    def __init__(self, name, mass, spin, parity):
        self.name = name
        self.mass = mass
        self.spin = spin
        self.parity = parity

    def spin_parity(self):
        return "${}^{}$".format(self.spin, self.parity)

    def __str__(self):
        return "{} ${}$ ({:.0f} MeV)".format(self.spin_parity(), self.name, self.mass)

class Resonance(Particle):
    def __init__(self, name, mass, spin, parity, color):
        super().__init__(name, mass, spin, parity)
        self.color = color

def multiply_parities(parity1, parity2):
    return '+' if parity1 == parity2 else '-'


def compute_allowed_spins(spin1, spin2):
    min_spin = abs(spin1 - spin2)
    max_spin = spin1 + spin2
    return list(range(min_spin, max_spin + 1))

MASS_DIFF = 50
def generate_latex_table(particles1, particles2, resonances=[]):
    is_symmetric = particles1 == particles2
    header = " & ".join([""] + [str(p) for p in particles2]) + " \\\\ \\hline"
    rows = [header]

    for i, p1 in enumerate(particles1):
        row = [str(p1)]
        for j, p2 in enumerate(particles2):
            if is_symmetric and j < i:
                row_content = "---"  
            else:
                combined_mass = p1.mass + p2.mass
                combined_parity = multiply_parities(p1.parity, p2.parity)
                spins = compute_allowed_spins(p1.spin, p2.spin)
                if len(spins) == 1:
                    spin_representation = "${}^{}$".format(spins[0], combined_parity)
                else:
                    spin_representation = "$({})^{}$".format(','.join(map(str, spins)), combined_parity)

                # Check if cell matches resonance criteria
                found = False
                cell_content = "{} ({:.0f} MeV)".format(spin_representation, combined_mass)

                for resonance in resonances:
                    if (resonance.mass - MASS_DIFF <= combined_mass <= resonance.mass + MASS_DIFF
                            and combined_parity == resonance.parity
                            and resonance.spin in spins):
                        
                        if not found:
                            # First matching resonance; color the entire content
                            cell_content = "\\textcolor{{{}}}{{{}}}".format(resonance.color, cell_content)
                            found = True
                        else:
                            # Subsequent matching resonances; only color the specific spin value
                            colored_spin = "\\textcolor{{{}}}{{{}}}".format(resonance.color, resonance.spin)
                            cell_content = cell_content.replace(str(resonance.spin), colored_spin)

                row_content = cell_content
            row.append(row_content)
        rows.append(" & ".join(row) + " \\\\ \\hline" )

    table = "\\begin{tabular}{" + "|c" * (len(particles2) + 1) + "|}\n"
    table += " \\hline \\hline \n"
    table += "\n".join(rows)
    table += " \\hline"
    table += "\n\\end{tabular}"

    return table

# From PDG
particlesD0 = [Particle("D^0", 1864.84, 0, "-"), 
               Particle("D^*", 2006.85, 1, "-"),  
               Particle("D_0^*(2300)^0", 2343, 0, "+"),  
               Particle("D_{1}(2420)^0", 2422.1, 1, "+"), 
               Particle("D_{1}(2430)^0", 2412, 1, "+"),   
               Particle("D^*_{2}(2460)^0", 2461, 2, "+"), 
               Particle("D_{0}(2550)^0", 2549, 0, "-"),  #
               Particle("D^*_{1}(2600)^0", 2672, 1, "-"), # 
               Particle("D^*_{2}(2740)^0", 2747, 2, "-"), 
               Particle("D^*_{3}(2750)^0", 2763.1, 3, "-"),  
               Particle("D^*_{1}(2760)^0", 2781, 1, "-")  #              
               ]
particlesD = [Particle("D^\pm", 1864.66, 0, "-"), 
               Particle("D^{*\pm}", 2010.26, 1, "-"),  
               Particle("D_0^*(2300)^\pm", 2343, 0, "+"),  
               Particle("D_{1}(2420)^\pm", 2422.1, 1, "+"), 
               Particle("D_{1}(2430)^\pm", 2412, 1, "+"),   
               Particle("D^*_{2}(2460)^\pm", 2461, 2, "+"), 
               #Particle("D^{*\pm}", 2637, , ""), 
               Particle("D^*_{3}(2750)^\pm", 2763.1, 3, "-"), 
               ]
particlesDs = [Particle("D_{s}^\pm", 1968.35, 0, "-"), 
               Particle("D_{s}^{*\pm}", 2112.2, 1, "-"),  
               Particle("D_{s0}^*(2317)^\pm", 2317.8, 0, "+"),  
               Particle("D_{s1}^*(2460)^\pm", 2459.5, 1, "+"),  
               Particle("D_{s1}^*(2536)^\pm", 2535.11, 1, "+"),  
               Particle("D_{s2}^*(2573)^\pm", 2569.1, 2, "+"),  
               Particle("D_{s0}^*(2590)^\pm", 2591, 0, "-"), #             
               Particle("D_{s1}^*(2700)^\pm", 2714, 1, "-"),  
               Particle("D_{s1}^*(2860)^\pm", 2859, 1, "-"), #
               Particle("D_{s3}^*(2860)^\pm", 2860, 3, "-"),  
               ]

particlesPsi = [
    Particle("\eta_c(1S)", 2983.9, 1, "-"),
    Particle("J/\psi(1S)", 3096.9, 1, "-"),
    Particle("\chi_{c0}(1P)", 3414.8, 1, "-"),
    Particle("\chi_{c1}(1P)", 3510.7, 1, "-"),
    Particle("\chi_{c2}(1P)", 3556.2, 1, "-"),
    Particle("\eta_c(2S)", 3639.0, 1, "-"),
    Particle("\psi(2S)", 3686.0, 1, "-"),
    Particle("\psi(3770)", 3770.0, 1, "-"),
    Particle("\psi(4040)", 4040.0, 1, "-"),
    Particle("\psi(4160)", 4160.0, 1, "-"),
    Particle("\psi(4415)", 4415.0, 1, "-")
]

particlesOther = [
    Particle("K^+", 493.7, 0, "-"),  
    Particle("K^*(892)^0", 895.5, 1, "-"),  
    Particle("\\pi^+", 139.6, 0, "-"),  
    Particle("\\rho(770)^0", 775.5, 1, "-")  
]

# From amp fit
resonancesX = [Resonance("\psi(4360)", 4372, 1, "-", "red"), Resonance("X(S)", 4482, 0, "+", "blue"), Resonance("X(A)", 4688, 1, "+", "green"), Resonance("X^2(S)", 4673, 0, "+", "orange"), Resonance("X(V)", 4778, 1, "-", "cyan")]
resonancesZ = [Resonance("X(4055)", 4054, 1, "-", "red"), Resonance("Z(A)", 4278, 1, "+", "blue"), Resonance("Z^2(A)", 4481, 1, "+", "green")]
resonancesXs = [Resonance("X_{s}(S)", 4410, 0, "+", "red"), Resonance("X_{s}(A)", 4576, 1, "+", "blue"), Resonance("X^2_{s}(A)", 4935, 1, "+", "green"), Resonance("X_{s}(V)", 5097, 1, "-", "orange")]
resonancesZs = [Resonance("Z_{s}(4000)", 4003, 1, "+", "red")]

# Create tables
# latex_table = generate_latex_table(particlesDs, particlesD, resonances=resonancesXs)
# print("\n" + latex_table + "\n")
# resonance_strs = ["\\textcolor{" + resonance.color + "}{$" + resonance.name + "$ (" + str(resonance.mass) + " MeV)}" for resonance in resonancesXs]
# print("Thresholds with the appropriate $J^P$ and within {} MeV of {} are highlighted.".format(MASS_DIFF, ', '.join(resonance_strs)))

# latex_table = generate_latex_table(particlesD0, particlesD0, resonances=resonancesX)
# print("\n" + latex_table + "\n")
# resonance_strs = ["\\textcolor{" + resonance.color + "}{$" + resonance.name + "$ (" + str(resonance.mass) + " MeV)}" for resonance in resonancesX]
# print("Thresholds with the appropriate $J^P$ and within {} MeV of {} are highlighted.".format(MASS_DIFF, ', '.join(resonance_strs)))

# latex_table = generate_latex_table(particlesD, particlesD, resonances=resonancesX)
# print("\n" + latex_table + "\n")
# resonance_strs = ["\\textcolor{" + resonance.color + "}{$" + resonance.name + "$ (" + str(resonance.mass) + " MeV)}" for resonance in resonancesX]
# print("Thresholds with the appropriate $J^P$ and within {} MeV of {} are highlighted.".format(MASS_DIFF, ', '.join(resonance_strs)))

# latex_table = generate_latex_table(particlesDs, particlesDs, resonances=resonancesX)
# print("\n" + latex_table + "\n")
# resonance_strs = ["\\textcolor{" + resonance.color + "}{$" + resonance.name + "$ (" + str(resonance.mass) + " MeV)}" for resonance in resonancesX]
# print("Thresholds with the appropriate $J^P$ and within {} MeV of {} are highlighted.".format(MASS_DIFF, ', '.join(resonance_strs)))

# latex_table = generate_latex_table(particlesDs, particlesD0, resonances=resonancesZs)
# print("\n" + latex_table + "\n")
# resonance_strs = ["\\textcolor{" + resonance.color + "}{$" + resonance.name + "$ (" + str(resonance.mass) + " MeV)}" for resonance in resonancesZs]
# print("Thresholds with the appropriate $J^P$ and within {} MeV of {} are highlighted.".format(MASS_DIFF, ', '.join(resonance_strs)))

# latex_table = generate_latex_table(particlesD, particlesD0, resonances=resonancesZs)
# print("\n" + latex_table + "\n")
# resonance_strs = ["\\textcolor{" + resonance.color + "}{$" + resonance.name + "$ (" + str(resonance.mass) + " MeV)}" for resonance in resonancesZ]
# print("Thresholds with the appropriate $J^P$ and within {} MeV of {} are highlighted.".format(MASS_DIFF, ', '.join(resonance_strs)))

# latex_table = generate_latex_table(particlesPsi, particlesOther, resonances=[])
# print("\n" + latex_table + "\n")


ffs = [0.383, 0.0157 ,0.232, 0.0045 ] # Belle fit 1
#ffs = [0.430, 0.0184 ,0.168, 0.00758]  # Belle fit 2
total = sum(ffs)
normalized_ffs = [x/total * 100 for x in ffs]
print(normalized_ffs)

bfs = [ffs[0], ffs[1] / 0.93 * 3./4., ffs[2] * 3./4., ffs[3] / 0.0153 * 1./3. ]
total = sum(bfs)
normalized_bfs = [x/total * 100 for x in bfs]
print(normalized_bfs)


ffs = [0.5707, 0.0889 ,0.136 + 0.0589, 0.00758]  # Amp fit
bfs = [ffs[0], ffs[1] * 3./4., ffs[2] * 3./4., 0.11, 0.03]
total = sum(bfs)
normalized_bfs = [x/total * 100 for x in bfs]
print(normalized_bfs)