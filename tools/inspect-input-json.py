import sys
import json
import gufe
import openfe


def write_sdf(smc):
    sdf = smc.to_sdf()
    with open(f"{smc.name}.sdf", 'w') as f:
        f.write(sdf)

with open(sys.argv[1]) as stream:
    r = json.load(stream, cls=gufe.tokenization.JSON_HANDLER.decoder)


# Get inputs
systemA = gufe.Component.from_dict(r['stateA'])
systemB = gufe.Component.from_dict(r['stateB'])
protocol = gufe.protocols.Protocol.from_dict(r['protocol'])

if r['mapping'] is not None:
    mapping = gufe.LigandAtomMapping.from_dict(r['mapping'])
else:
    mapping = None

print(systemA)
print(systemB)
print(protocol)

# Write out components and mapping
for comp in systemA.components.values():
    if isinstance(comp, gufe.ProteinComponent):
        comp.to_pdb_file('protein.pdb')

# Print out the ligands
for comp in systemA.components.values():
    if isinstance(comp, gufe.SmallMoleculeComponent):
        write_sdf(comp)

for comp in systemB.components.values():
    if comp not in systemA.components.values():
        write_sdf(comp)

if mapping is not None:
    mapping.draw_to_file('mapping.png')
