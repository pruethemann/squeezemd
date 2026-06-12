import pandas as pd


def main():
    overview_df = pd.read_parquet("/home/peter/caracara/Squeeze/md_overview.parquet")

    print(overview_df)


if __name__ == "__main__":
    main()
