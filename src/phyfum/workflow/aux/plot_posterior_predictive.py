import numpy as np
import pandas as pd
from scipy import linalg, stats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import click


def generate_state_var(S):
    return [(k, m) for m in range(S + 1) for k in range(S + 1) if k + m <= S]


def generate_rate_matrix(S, lam, mu, gamma, stateVar=None):
    if stateVar is None:
        stateVar = generate_state_var(S)

    RateMatrix = np.zeros((len(stateVar), len(stateVar)))

    for down, (k_down, m_down) in enumerate(stateVar):
        for across, (k, m) in enumerate(stateVar):
            if k == k_down - 1 and m == m_down:
                RateMatrix[down, across] = (S - m - k) * (k * lam / (S - 1) + 2 * mu)
            elif k == k_down and m == m_down - 1:
                RateMatrix[down, across] = m * (S - m - k) * lam / (S - 1)
            elif k == k_down + 1 and m == m_down - 1:
                RateMatrix[down, across] = k * (m * lam / (S - 1) + mu)
            elif k == k_down + 1 and m == m_down:
                RateMatrix[down, across] = k * ((S - m - k) * lam / (S - 1) + gamma)
            elif k == k_down and m == m_down + 1:
                RateMatrix[down, across] = m * (S - m - k) * lam / (S - 1)
            elif k == k_down - 1 and m == m_down + 1:
                RateMatrix[down, across] = m * (k * lam / (S - 1) + 2 * gamma)
            elif k == k_down and m == m_down:
                RateMatrix[down, across] = -(2 * ((k + m) * (S - m - k) + k * m) * lam / (S - 1) + (k + 2 * m) * gamma + (2 * S - (k + 2 * m)) * mu)

    return RateMatrix


def findProbDist(RateMatrix, InitialConditions, stateVar, age, S):

    ProbStates = linalg.expm(RateMatrix * age) @ InitialConditions
    ProbDist = np.zeros(2 * S + 1)
    for index, (k, m) in enumerate(stateVar):
        ProbDist[k + 2 * m] += ProbStates[index]

    return ProbDist


def runModel(S, lam, mu, gamma, age):

    stateVar = generate_state_var(S)
    RateMatrix = generate_rate_matrix(S, lam, mu, gamma, stateVar)
    InitialConditions = np.zeros(int(0.5 * (S + 1) * (S + 2)))
    InitialConditions[[0, -1]] = 0.5
    ProbDist = findProbDist(RateMatrix, InitialConditions, stateVar, age, S)

    return ProbDist


def beta_convert_params(mu, kappa):
    """
    Convert mean/dispersion parameterization of a beta distribution to the ones scipy supports

    """

    if np.any(kappa <= 0):
        raise Exception("kappa must be greater than 0")
    elif np.any(mu <= 0) or np.any(mu >= 1):
        raise Exception("mu must be between 0 and 1")

    alpha = kappa * mu
    beta = kappa * (1 - mu)

    return alpha, beta


def rescale_beta(beta, delta, eta):
    # Linear transform of beta values from between
    # 0 and 1 to between delta and eta
    return (eta - delta) * beta + delta


def add_noise(beta, delta, eta, kappa):
    beta_rescale = rescale_beta(beta, delta, eta)
    a, b = beta_convert_params(beta_rescale, kappa)

    return stats.beta.rvs(a, b)


def simulate_beta(params, S, age, N):
    lam, mu, gamma, delta, eta, kappa = params

    ProbDist = runModel(S, lam, mu, gamma, age)

    k_sample = np.random.choice(np.arange(0, 2 * S + 1), size=N, p=ProbDist)
    beta_sample = k_sample / (2 * S)

    return add_noise(beta_sample, delta, eta, kappa)


@click.command()
@click.option("--betas", help="Original csv file with the beta files", required=True, type=click.Path(exists=True))
@click.option("--log", help="BEAST log output", required=True, type=click.Path(exists=True))
@click.option("--samplesheet", help="Metadata CSV containing sample and age columns", required=True, type=click.Path(exists=True))
@click.option('--output', type=click.Path(exists=False), required=True, help='Output name (pdf)')
@click.option("--sample_col", help="Number of stem cells used in the model", required=True, type=str)
@click.option("--age_col", help="Number of stem cells used in the model", required=True, type=str)
def main(betas, log, samplesheet, output, sample_col, age_col):

    assert samplesheet.endswith(".csv") and betas.endswith(".csv") and log.endswith(".log"), "Unexpected input format. The betas and samplesheet should be a .csv file, while log should be a .log file"
    df = pd.read_csv(betas, index_col=0)
    df.columns = df.columns.astype(str)

    posterior = pd.read_csv(log, index_col=0, sep="\t", comment="#")
    # Get the column alignment.stemCells as a list of integer
    stemCells = posterior["alignment.stemCells"].values.astype(int)

    posterior = posterior[["flipflop.lambda", "flipflop.mu", "flipflop.gamma", "errorModel.deltaOffset", "errorModel.etaOffset", "errorModel.kappaScale"]]

    samplesheet = pd.read_csv(samplesheet)
    ages = samplesheet.groupby(sample_col)[age_col].max().to_dict()

    with PdfPages(output) as export_pdf:
        for sample in df.columns:
            y = df[sample]
            N = y.shape[0]
            P = posterior.shape[0]
            age = ages[sample]
            
            y_hat = np.empty((P, N))

            for p in range(P):
                params = posterior.iloc[p].values
                S = stemCells[p]
                y_hat[p, :] = simulate_beta(params, S, age, N)

            fig, ax = plt.subplots()
            plt.hist(y, np.linspace(0, 1, 101), density=True, alpha=0.4, linewidth=0)
            plt.hist(np.ravel(y_hat), np.linspace(0, 1, 101), density=True, alpha=0.4, linewidth=0)
            plt.legend(("Data", "Posterior Predictive"))
            plt.xlabel("Fraction Methylated (Beta)")
            plt.ylabel("Probability Density")
            plt.title(f"Sample ID: {sample}")
            sns.despine()
            plt.tight_layout()
            export_pdf.savefig(dpi=600)
            plt.close()


if __name__ == "__main__":
    main()
