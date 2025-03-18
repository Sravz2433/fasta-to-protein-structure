
from modeller import *
from modeller.automodel import *  # Load 'automodel' for homology modeling

def run_modeller(fasta_sequence, template_pdb, alignment_file="alignment.ali", output_model="model.pdb"):
    """
    Generate a homology model using MODELLER.
    
    Args:
        fasta_sequence (str): Protein sequence in FASTA format.
        template_pdb (str): PDB file of the homologous template.
        alignment_file (str): Alignment file (.ali) used for MODELLER.
        output_model (str): Output filename for the predicted model.
    """

    # Step 1: Create MODELLER environment
    env = Environ()
    env.io.atom_files_directory = ['.', './templates']

    # Step 2: Define homology modeling class
    class MyModel(AutoModel):
        def special_patches(self, aln):
            self.rename_segments(segment_ids=['A'], renumber_residues=[1])

    # Step 3: Run Homology Modeling
    a = MyModel(env, 
                alnfile=alignment_file,       # Alignment file
                knowns=template_pdb,          # Template PDB ID
                sequence="target",            # Target sequence name
                assess_methods=(assess.DOPE, assess.GA341))

    a.starting_model = 1
    a.ending_model = 5  # Generate 5 models and choose the best one

    a.make()

    # Step 4: Select the best model
    ok_models = [x for x in a.outputs if x['failure'] is None]
    best_model = min(ok_models, key=lambda x: x['DOPE score'])
    best_model_name = best_model['name']

    print(f"✅ Best model selected: {best_model_name}")

    # Rename and save the final model
    import shutil
    shutil.move(best_model_name, output_model)
    print(f"✅ Model saved as: {output_model}")

# Example: Run the function
run_modeller("MKTLVVVTIVCLDLGYTPE", "template.pdb")
