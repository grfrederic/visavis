from visavis import VisAVisClient
from parameters.default import parameters


# initialize a new client
client = VisAVisClient(
    infection_bin='../target/release/vis-a-vis',  # path to compiled binary
    sim_root='outdir',  # will be created if it doesn't exist
)

result = client.run(
    parameters_json=parameters,
    protocol_file_path='../protocols/default.protocol',
    clean_up=False,  # don't remove files from outdir/ after simulation
    images=True,     # save images at specified time points...
    annotate=True,   # ...and annotate them
)

# access simulation result easily
print("\nstates and neighbors as pandas dataframes:")
print(result.states)
print(result.neighbors)
