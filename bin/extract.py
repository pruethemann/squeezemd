import pandas as pd

def parse_lipophilic(parts):
        interaction_info = parts[0].split()
        donor_acceptor = parts[1].strip().split()

        receptor_resname = donor_acceptor[1]
        receptor_resid = donor_acceptor[2]

        ligand_resname = donor_acceptor[-2]
        ligand_resid = donor_acceptor[-1]
        distance = float(interaction_info[2].split("=")[1])
        energy = float(interaction_info[3].split("=")[1])

        interaction = {
            "Interaction Type": 'lipophilic',
            "Distance (r)": distance,
            "Energy (e)": energy,
            'receptor_resname' : receptor_resname,
            'receptor_resid' : receptor_resid,
            'ligand_resname' : ligand_resname,
            'ligand_resid' : ligand_resid,
        }

        return interaction



def parse_hbonds(parts):

        interaction_info = parts[0].split()
        donor_acceptor = parts[1].strip().split()

        receptor_resname = donor_acceptor[1]
        receptor_resid = donor_acceptor[2]

        ligand_resname = donor_acceptor[-2]
        ligand_resid = donor_acceptor[-1]

        distance = float(interaction_info[2].split("=")[1])
        angle = float(interaction_info[3].split("=")[1])
        energy = float(interaction_info[4].split("=")[1])


        marked = "marked as salt-bridge" in line

        interaction = {
            "Interaction Type": 'H-bond',
            "Distance (r)": distance,
            "Angle (a)": angle,
            "Energy (e)": energy,
            'receptor_resname' : receptor_resname,
            'receptor_resid' : receptor_resid,
            'ligand_resname' : ligand_resname,
            'ligand_resid' : ligand_resid,
            "Marked as Salt-Bridge": marked
        }

        return interaction

# Parse the input data into a pandas DataFrame
def parse():
    data = []

    with open("output.txt", 'r') as file:
        for line in file:

            if line.startswith("Lipo_EXT:"):
                parts = line.split("  !  ")
                interaction = parse_lipophilic(parts)
                data.append(interaction)

            if line.startswith("HB_EXT:"):
                parts = line.split("  !  ")
                interaction = parse_lipophilic(parts)
                data.append(interaction)

    return pd.DataFrame(data)

# Convert the input data into a DataFrame
interactions = parse()

# Save the DataFrame to a CSV file
interactions.to_csv("interactions.csv", index=False)

print("Data saved")
