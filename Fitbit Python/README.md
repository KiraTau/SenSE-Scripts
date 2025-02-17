# Fitbit Data Download

This folder contains an example script that allows a user to download their own data. It serves only as a utility for easily downloading the data for research purposes and helps walk through the basic concepts needed to build a fully-fledged app that accesses Fitbit data. [See the SenSE App here](https://github.com/Chukwuemeka-Ike/SenSEApp).

# Attribution
1. The main data_download script is based on the work by Michael Galarnyk in this [Medium article](https://towardsdatascience.com/using-the-fitbit-web-api-with-python-f29f119621ea).
2. The contents of the python_fitbit folder are copied from the [python-fitbit project](https://github.com/orcasgit/python-fitbit). Slight modifications have been made to allow one use Python 3.8.10 since the original module is no longer updated.

# Usage
All testing was done with Python 3.8.10, so I can only guarantee the script will function with that version. Please modify the commands below 

1. Install required modules.
```bash
pip install -r requirements.txt
pip install -r requirements/base.txt
```
2. Update the Client ID and Secret in data_download.py.
```python
# **************************** MODIFY THIS. ****************************
# Set the Client ID + Secret to enable access to the user's data.
CLIENT_ID = ""
CLIENT_SECRET = ""
```
3. Set the strings for the time series of interest. I've included the values for heart rate and steps. More information can be found on Fitbit's API documentation ([helpful starting point](https://dev.fitbit.com/build/reference/web-api/heartrate-timeseries/get-heartrate-timeseries-by-date-range/)).
4. Set the start and end dates of interest. Note that some Fitbit time series do not allow date ranges. Heart rate and steps do.
4. Run the script
```bash
python data_download.py
```
If successful, the script plots the time series and saves the results to a CSV in the same directory it was run from.

