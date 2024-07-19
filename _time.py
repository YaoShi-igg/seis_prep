from obspy import UTCDateTime

def utc2str(utc_time):
    """
    Convert UTCDateTime to string format.
    Args:
        - utc_time (UTCDateTime): The UTCDateTime object to convert.
    Returns:
        - str: The time as a string in the format yyyymmddhhmmssmmm.
    """
    # Format date as yyyymmdd
    date = ''.join(str(utc_time).split('T')[0].split('-'))
    # Format time as hhmmssmmm
    time = ''.join(str(utc_time).split('T')[1].split(':'))[:10]

    time_all = date + time
    str_time = ''.join(time_all.split('.'))

    return str_time

def utc2int(utc_time):
    """
    Convert UTCDateTime to integer format.
    Args:
        - utc_time (UTCDateTime): The UTCDateTime object to convert.
    Returns:
        - int: The time as an integer in the format yyyymmddhhmmssmmm.
    """
    int_time = int(utc2str(utc_time))
    return int_time

def str2utc(str_time):
    """
    Convert a string time back to UTCDateTime.
    Args:
        - str_time (str): Time as a string in the format yyyymmddhhmmssmmm.
    Returns:
        - UTCDateTime: The UTCDateTime object representing the given time.
    """
    year = int(str_time[0:4])
    month = int(str_time[4:6])
    day = int(str_time[6:8])
    hour = int(str_time[8:10])
    minute = int(str_time[10:12])
    second = int(str_time[12:14])
    microsecond = int(str_time[14:17]) * 1000  # Convert milliseconds to microseconds

    utc_time = UTCDateTime(year, month, day, hour, minute, second, microsecond)
    return utc_time

def add_date(date, n):
    """
    Increment a given date by n days.
    Args:
        - date (str): Date in the format YYYYMMDD.
        - n (int): Number of days to add.
    Returns:
        - str: New date in the format YYYYMMDD.
    """
    # Extract year, month, and day from the input date
    year = int(date[:4])
    month = int(date[4:6])
    day = int(date[6:8])
    
    # Define the number of days in each month
    per_month_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    
    for _ in range(n):
        # Check for leap year
        if (year % 400 == 0) or (year % 4 == 0 and year % 100 != 0):
            per_month_days[1] = 29
        else:
            per_month_days[1] = 28
        
        # Increment day
        if day < per_month_days[month - 1]:
            day += 1
        else:
            day = 1
            if month < 12:
                month += 1
            else:
                month = 1
                year += 1
    
    # Return the new date in YYYYMMDD format
    return f'{year:04d}{month:02d}{day:02d}'    

def get_next_start_of_day(utc_time):
    """
    Given a UTCDateTime object, returns the start of the next day (00:00:00)
        if the time part is not already at the start of the day.
    Args:
        - utc_time (UTCDateTime): The input UTCDateTime object.

    Returns:
        - UTCDateTime: A new UTCDateTime object set to the start of the next day.
    """
    # Check if the time part is not already at the start of the day
    if utc_time.hour != 0 or utc_time.minute != 0 or utc_time.second != 0:
        # Move to the next day and set the time to 00:00:00
        next_day = utc_time + (24 * 3600)  # Add 24 hours in seconds
        next_day = UTCDateTime(year=next_day.year, month=next_day.month, day=next_day.day)
    else:
        next_day = utc_time
    
    return next_day


if __name__ == '__main__':

    utc_time = UTCDateTime(2021, 9, 17, 10, 10, 23, 13)
    str_time = utc2str(utc_time)
    int_time = utc2int(utc_time)
    print("Original UTC Time:", utc_time)
    print("Converted String:", str_time)
    print("Converted Integer:", int_time)
    print("Converted Back to UTC:", str2utc(str_time))

    initial_date = '20210529'  # YYYYMMDD format
    days_to_add = 10
    new_date = add_date(initial_date, days_to_add)
    print(f"New date after adding {days_to_add} days: {new_date}")
