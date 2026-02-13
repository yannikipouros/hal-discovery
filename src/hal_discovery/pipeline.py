from .steps import step1_family_space, step2_interpro_retrieval

def run(ctx: RunContext) -> None:
    # (same validation as before...)

    families_csv = step1_family_space.run(ctx)

    # Part 2
    step2_interpro_retrieval.run(ctx, families_csv)

    print("\n✅ Part 2 complete.")
    print("Next: Part 3 will compile the retrieved accession tables into one dataframe, deduplicate, and filter.")
