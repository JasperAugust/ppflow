from ppflow.datasets.peptide_experiment import PairDataset

# Create an instance of PairDataset
dataset = PairDataset(
    data_dir='./dataset/nanomed/',
    processed_dir='./dataset/nanomed/processed',
    split_seed=2024,
    reset=True  # Set this to True if you want to regenerate the splits
)

# Now you can use the dataset for testing and experiments
print(f"Total number of samples: {len(dataset)}")

# Access a sample
sample = dataset[0]
print(f"Sample keys: {sample.keys()}")