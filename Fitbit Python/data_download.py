import python_fitbit.gather_keys_oauth2 as Oauth2
from python_fitbit.fitbit.api import Fitbit
import pandas as pd     
import json
import matplotlib.pyplot as plt

# **************************** MODIFY THIS. ****************************
# Set the Client ID + Secret to enable access to the user's data.
CLIENT_ID = "23R4R2"
CLIENT_SECRET = "1aa607e465e9306c30b9ebfeebf5bbcb"

# Uncomment the data you want.
data = "activities/heart"
data_response = "activities-heart-intraday"
filename = "heart"

# data = "activities/steps"
# data_response = "activities-steps-intraday"
# filename = "steps"

# Set the dates for the data request.
start_date = pd.to_datetime("2024/04/04")
end_date = pd.to_datetime("2024/04/08")

# print()
# print(f"Start time: {start_date}")
# print(f"End time: {end_date}")
# print()

# Start the server with the client credentials.
server = Oauth2.OAuth2Server(CLIENT_ID, CLIENT_SECRET)
server.browser_authorize()
ACCESS_TOKEN = str(server.fitbit.client.session.token['access_token'])
REFRESH_TOKEN = str(server.fitbit.client.session.token['refresh_token'])
auth2_client = Fitbit(
        CLIENT_ID,
        CLIENT_SECRET,
        oauth2=True,
        access_token=ACCESS_TOKEN,
        refresh_token=REFRESH_TOKEN
)

#
date_list = []
df_list = []
all_dates = pd.date_range(start=start_date, end=end_date)

for one_day in all_dates:
    one_day = one_day.date().strftime("%Y-%m-%d")
    one_day_data = auth2_client.intraday_time_series(
        data,
        base_date=one_day,
        detail_level="1min"
    )

    one_dayDf = pd.DataFrame(one_day_data[data_response]["dataset"])
    date_list.append(one_day)
    df_list.append(one_dayDf)

# Clean it up.
final_df_list = []
for date, df in zip(date_list, df_list):
    if len(df) == 0:
        continue
    df.loc[:, "date"] = pd.to_datetime(date)
    # df.loc[:, "datetime"] = df.time
    df['datetime'] = df['date'] + pd.to_timedelta(df['time'])
    final_df_list.append(df)
print(final_df_list)
final_df = pd.concat(final_df_list, axis=0)

# Print the head for sanity check.
print(final_df.info())
print(final_df.head())

# Plot the values.
final_df.plot(
    kind="line",
    x="datetime",
    y="value"
)
plt.show()

# Save the data.
final_df.to_csv(filename + '.csv', index=False)