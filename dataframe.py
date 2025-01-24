"""
dataframe.py
============

This code provides functionalities for analyzing experiment data from an experimentation platform. 
It focuses on merging datasets, calculating metrics like QWRI (Quality-Weighted Revenue Impact), 
and performing statistical analyses to assess the effectiveness of experiments.

Key Features:
-------------
- **Data Integration**:
    - Merges multiple experiment-related datasets, including experiment details, tracking data, and user exposure data.
- **Metric Calculations**:
    - Computes key metrics such as:
        - **Gross Projected Future Impact**: Projects the full-scale impact of experiments based on audience fractions.
        - **Future Estimated Impact**: Adjusts for quality using a discount factor that accounts for statistical power and sample ratio mismatches.
        - **QWRI**: A comprehensive metric combining revenue impact and statistical success.
- **Statistical Analysis**:
    - Computes correlation matrices between variables to uncover trends.
    - Calculates the statistical success rate of experiments.

Dependencies:
-------------
- pandas: For data manipulation and analysis.    
"""

# IMPORT LIBRARIES
import pandas as pd

# CREATING PANDAS DATAFRAME FROM THE DATA TABLE
data = {
    "Month": ["Apr-24", "May-24", "Jun-24", "Jul-24", "Aug-24", "Sep-24", "Oct-24"],
    "Experiments Launched": [40, 74, 67, 69, 85, 92, 108],
    "Estimated Revenue Impact": [1197380, 24837, -42, 601995, 54100, 353, 214998],
    "Estimated Risk Mitigated": [32810, 30858, -45393, 0, -398536, -6470, -1203],
}

df = pd.DataFrame(data)
df["Experiments Launched"] = df["Experiments Launched"].astype(int)
df["Estimated Revenue Impact"] = df["Estimated Revenue Impact"].astype(int)
df["Estimated Risk Mitigated"] = df["Estimated Risk Mitigated"].astype(int)

# CREATING PANDAS DATAFRAME FROM THE EXPERIMENTS TABLE
experiments = {
    "Month": ["Apr-24", "May-24", "Jun-24", "Jul-24", "Aug-24", "Sep-24", "Oct-24"],
    "Completed Experiments": [68, 48, 55, 73, 54, 66, 93],
    "Winning Treatment": [36, 24, 25, 38, 29, 17, 40],
    "Winning Control": [32, 24, 30, 35, 25, 49, 53],
}

df2 = pd.DataFrame(experiments)
df2["Completed Experiments"] = df2["Completed Experiments"].astype(int)
df2["Winning Treatment"] = df2["Winning Treatment"].astype(int)
df2["Winning Control"] = df2["Winning Control"].astype(int)

combined_df = pd.merge(df, df2, on="Month")
# PRINT BOTH DATAFRAMES MERGED INTO ONE ON THE "MONTH" COLUMN
# SIMILAR TO SQL JOIN
print('COMBINED DATAFRAMES\n')
print(combined_df)
print()

# FUNCTION TO CALCULATE DISCOUNT FACTOR AS PROVIDED IN THE ASSIGNMENT:
def discount_factor(primary_metric_power: float, sample_ratio_mismatch: bool) -> float:
    """
    Calculate the discount factor based on primary metric power and sample ratio mismatch.

    The discount factor adjusts the weight of an experiment's impact based on its statistical 
    power and whether a sample ratio mismatch (SRM) exists. If SRM is present, the discount 
    factor is set to 0. Otherwise, the discount factor is calculated as:
    
        (2 * primary_metric_power) / (1 + primary_metric_power)

    Args:
        primary_metric_power (float): The statistical power of the primary metric, ranging 
                                      between 0 and 1. Higher values indicate better 
                                      experiment reliability.
        sample_ratio_mismatch (bool): A flag indicating whether there is a sample ratio 
                                       mismatch in the experiment. True means there is an 
                                       imbalance, and the discount factor is set to 0.

    Returns:
        float: The calculated discount factor, ranging from 0 to 1. A value of 0 indicates 
               either SRM is present or the power is insufficient for meaningful impact.

    Examples:
        >>> discount_factor(0.8, False)
        0.8888888888888888

        >>> discount_factor(0.8, True)
        0

        >>> discount_factor(0.5, False)
        0.6666666666666666
    """
    return 0 if sample_ratio_mismatch else (2 * primary_metric_power) / (1 + primary_metric_power)

# CALCULATE CORRELATIONS BETWEEN ALL VARIABLES AND PRINT THE MATRIX
correlation_matrix = combined_df[["Experiments Launched", "Estimated Revenue Impact", "Estimated Risk Mitigated",
                                  "Completed Experiments", "Winning Treatment", "Winning Control"]].corr()

print("CORRELATION MATRIX\n")
print(correlation_matrix)
print()

# DATA TABLES

"""
EXPERIMENT
Description:
  exp_id : unique identifier for each experiment
  start_date: start date of an experiment
  end_date: end date of an experiment
  description: test description
  fraction: fraction of users exposed to the test
"""

experiment_data = {
    "exp_id": [1, 2, 3],
    "start_date": ["2022-01-01", "2022-05-17", "2023-10-22"],
    "end_date": ["2022-02-10", "2022-06-21", "2023-11-15"],
    "description": ["Test new button X", "Test new onboarding flow", "Tet new feature B"],
    "fraction": [10, 16, 20],
}

# Creating the "experiment" DataFrame
experiment = pd.DataFrame(experiment_data)

"""
EXPERIMENT TRACKING
Description:
  timestamp: timestamp of the event
  user_id: unique identifier for each user
  exp_id: unique identifier for each experiment
  exp_variant: experiment variant
  event_type: events measured within the experiment, such as purchases,
  subscription renewals, unsubscribes, etc.
  value: numeric value associated with the event
"""

experiment_tracking_data = {
    "timestamp": ["2022-01-02 15:30:00", "2022-01-03 18:39:00", "2022-05-03 07:01:18", "2023-10-23 19:14:45"],
    "user_id": ["user1", "user2", "user3", "user54"],
    "exp_id": [1, 1, 1, 3],
    "exp_variant": ["a", "b", "b", "a"],
    "event_type": ["purchase", "purchase", "unsubscribe", "purchase"],
    "value": [20.00, 30.00, 1.00, 10.00],
}

# Creating the "experiment_tracking" DataFrame
experiment_tracking = pd.DataFrame(experiment_tracking_data)


"""
EXPERIMENT USERS
Description:
Table contains the list of user ids that were exposed to an experiment.
  user_id: unique identifier for each user
  exp_id: experiment id that the user was exposed to
  variant_id: experiment variant that the user was assigned to
"""

experiment_users_data = {
    "user_id": ["user1", "user2", "user8"],
    "exp_id": [1, 1, 4],
    "variant": ["a", "b", "c"],
}

# Creating the "experiment_users" DataFrame
experiment_users = pd.DataFrame(experiment_users_data)


# CALCULATIONS

# FUNCTION TO CALCULATE GROSS PROJECTED FUTURE IMPACT AS PER THE ASSIGNMENT
def calculate_gross_projected_future_impact(exp_id, experiment, experiment_tracking, experiment_users):
    """
    Calculate the Gross Projected Future Impact of a single experiment.

    Args:
        exp_id (int): The experiment ID to calculate the impact for.
        experiment (pd.DataFrame): DataFrame containing experiment details (start_date, end_date, fraction).
        experiment_tracking (pd.DataFrame): DataFrame containing tracking data (event_type, value).
        experiment_users (pd.DataFrame): DataFrame containing user exposure data for experiments.

    Returns:
        float: Gross Projected Future Impact of the experiment.
    """
    # Get the experiment details
    exp_details = experiment[experiment["exp_id"] == exp_id]
    if exp_details.empty:
        raise ValueError(f"Experiment ID {exp_id} not found in the experiment DataFrame.")
   
    start_date = pd.to_datetime(exp_details["start_date"].values[0])
    end_date = pd.to_datetime(exp_details["end_date"].values[0])
    audience_fraction = exp_details["fraction"].values[0]
   
    # Calculate duration of the experiment
    duration_of_impact = (end_date - start_date).days + 1  # Include the start and end day
   
    # Filter tracking data for the experiment and calculate daily impact per unit
    tracking_data = experiment_tracking[experiment_tracking["exp_id"] == exp_id]
    if tracking_data.empty:
        raise ValueError(f"No tracking data found for Experiment ID {exp_id}.")
    
    total_event_value = tracking_data["value"].sum()  # Sum the values for all events
    daily_impact_per_unit = total_event_value / duration_of_impact
    
    # Calculate total exposed units from experiment_users
    total_exposed_units = experiment_users[experiment_users["exp_id"] == exp_id]["user_id"].nunique()
    if total_exposed_units == 0:
        raise ValueError(f"No users exposed to Experiment ID {exp_id}.")
    
    # Calculate Gross Projected Future Impact
    gross_projected_future_impact = (
        daily_impact_per_unit * total_exposed_units * duration_of_impact / audience_fraction
    )
    
    return gross_projected_future_impact

# Example usage
result = calculate_gross_projected_future_impact(1, experiment, experiment_tracking, experiment_users)
print(f"Gross Projected Future Impact for Experiment ID 1: {result}")


# FUNCTION TO CALCULATE FUTURE ESTIMATED IMPACT AS PER THE ASSIGNMENT
def calculate_future_estimated_impact(exp_id, experiment, experiment_tracking, experiment_users, primary_metric_power, sample_ratio_mismatch):
    """
    Calculate the Future Estimated Impact of a single experiment.

    Args:
        exp_id (int): The experiment ID to calculate the impact for.
        experiment (pd.DataFrame): DataFrame containing experiment details (start_date, end_date, fraction).
        experiment_tracking (pd.DataFrame): DataFrame containing tracking data (event_type, value).
        experiment_users (pd.DataFrame): DataFrame containing user exposure data for experiments.
        primary_metric_power (float): The statistical power of the primary metric (0 to 1).
        sample_ratio_mismatch (bool): Whether there is a sample ratio mismatch (True or False).

    Returns:
        float: Future Estimated Impact of the experiment.
    """
    # Calculate Gross Projected Future Impact
    exp_details = experiment[experiment["exp_id"] == exp_id]
    if exp_details.empty:
        raise ValueError(f"Experiment ID {exp_id} not found in the experiment DataFrame.")
    
    start_date = pd.to_datetime(exp_details["start_date"].values[0])
    end_date = pd.to_datetime(exp_details["end_date"].values[0])
    audience_fraction = exp_details["fraction"].values[0]
    duration_of_impact = (end_date - start_date).days + 1  # Include the start and end day
    
    tracking_data = experiment_tracking[experiment_tracking["exp_id"] == exp_id]
    if tracking_data.empty:
        raise ValueError(f"No tracking data found for Experiment ID {exp_id}.")
    
    total_event_value = tracking_data["value"].sum()
    daily_impact_per_unit = total_event_value / duration_of_impact
    
    total_exposed_units = experiment_users[experiment_users["exp_id"] == exp_id]["user_id"].nunique()
    if total_exposed_units == 0:
        raise ValueError(f"No users exposed to Experiment ID {exp_id}.")
    
    gross_projected_future_impact = (
        daily_impact_per_unit * total_exposed_units * duration_of_impact / audience_fraction
    )
    
    # Calculate Discount Factor
    discount = discount_factor(primary_metric_power, sample_ratio_mismatch)
    
    # Calculate Future Estimated Impact
    future_estimated_impact = gross_projected_future_impact * discount
    
    return future_estimated_impact

# Example usage
# Assuming experiment, experiment_tracking, and experiment_users are predefined DataFrames
primary_metric_power = 0.8
sample_ratio_mismatch = False

result = calculate_future_estimated_impact(1, experiment, experiment_tracking, experiment_users, primary_metric_power, sample_ratio_mismatch)
print(f"Future Estimated Impact for Experiment ID 1: {result}")


# CALCULATE A NEW METRIC (QWRI)
def calculate_qwri(experiment, experiment_tracking, experiment_users, primary_metric_power, sample_ratio_mismatch):
    """

    Calculate the Quality-Weighted Revenue Impact (QWRI) for the platform. (Higher = Better)
    QWRI = (total revenue impact from successful experiments x statistical success rate) / 100
    Statistical Success Rate = (number of statistically significant experiments / total experiments launched) x 100
    Statistically significant experiments are determined by their p-value.

    Args:
        experiment (pd.DataFrame): DataFrame containing experiment details.
        experiment_tracking (pd.DataFrame): DataFrame containing tracking data.
        experiment_users (pd.DataFrame): DataFrame containing user exposure data.
        primary_metric_power (float): The statistical power of the primary metric (0 to 1).
        sample_ratio_mismatch (bool): Whether there is a sample ratio mismatch (True or False).

    Returns:
        float: Quality-Weighted Revenue Impact (QWRI).
    """
    total_revenue_impact = 0
    successful_experiments = 0

    # Iterate through each experiment
    for exp_id in experiment["exp_id"]:
        try:
            # Calculate Gross Projected Future Impact
            gross_impact = calculate_gross_projected_future_impact(
                exp_id, experiment, experiment_tracking, experiment_users
            )
            discount = discount_factor(primary_metric_power, sample_ratio_mismatch)
            future_estimated_impact = gross_impact * discount

            # Increment successful experiments count
            successful_experiments += 1

            # Add to total revenue impact
            total_revenue_impact += future_estimated_impact

        except ValueError as e:
            # Log skipped experiments
            print(f"Skipping Experiment ID {exp_id} due to error: {e}")
            continue

    # Calculate Statistical Success Rate
    total_experiments_launched = len(experiment["exp_id"])
    statistical_success_rate = (successful_experiments / total_experiments_launched) * 100

    # Calculate QWRI
    if total_revenue_impact == 0:
        return 0  # Avoid division by zero
    qwri = (total_revenue_impact * statistical_success_rate) / 100
    return qwri

# TEST THE FUNCTION:

# Example inputs for primary_metric_power and sample_ratio_mismatch
primary_metric_power = 0.8  # Assumed statistical power (adjust as needed)
sample_ratio_mismatch = False  # Assume no sample ratio mismatch for simplicity

# Call the function and calculate QWRI
qwri = calculate_qwri(
    experiment=experiment, 
    experiment_tracking=experiment_tracking, 
    experiment_users=experiment_users, 
    primary_metric_power=primary_metric_power, 
    sample_ratio_mismatch=sample_ratio_mismatch
)

# Print the result
print(f"Quality-Weighted Revenue Impact (QWRI): {qwri:.2f}")