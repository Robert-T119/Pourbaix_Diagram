# Pourbaix_Diagram

## Overview

This is a Django-based web application that provides various functionalities related to Pourbaix Diagrams. It includes multiple apps such as `copper`, `nickelammonia`, `nickelca`, `nickelcitrate`, and `nickelpure`.

## Installation

### Prerequisites

- Python 3.x
- pip
- Virtual environment (optional but recommended)

### Steps

1. **Clone the Repository**
    ```bash
    git clone https://github.com/Robert-T119/Pourbaix_Diagram.git
    ```

2. **Navigate to the Project Directory**
    ```bash
    cd Pourbaix_Diagram
    ```

3. **Create a Virtual Environment (Optional)**
    ```bash
    python3 -m venv venv
    ```

4. **Activate the Virtual Environment (Optional)**
    - On Windows: `.\venv\Scripts\activate`
    - On macOS and Linux: `source venv/bin/activate`

5. **Install Dependencies**
    ```bash
    pip install -r requirements.txt
    ```

6. **Run Migrations**
    ```bash
    python manage.py migrate
    ```

7. **Run the Server**
    ```bash
    python manage.py runserver
    ```

## How to Run

After completing the installation, you can run the server using the following command:

```bash
python manage.py runserver
```

Then, open your web browser and navigate to [http://127.0.0.1:8000/](http://127.0.0.1:8000/).

## Features

- **Copper App**: Functions related to copper in Pourbaix Diagrams.
- **Nickelammonia App**: Manages nickel-ammonia interactions in Pourbaix Diagrams.
- **Nickelca App**: Focuses on nickel-calcium interactions.
- **Nickelcitrate App**: Handles nickel-citrate interactions in Pourbaix Diagrams.
- **Nickelpure App**: Features for pure nickel in Pourbaix Diagrams.


