Application that accepts a FASTQ format file of DNA sequence reads and assembles the reads into a single DNA sequence. This sequence is then saved as a FASTA format file for downloading. The application is designed to be deployed on Amazon Web Services as an Elastic Beanstalk application and was developed for Brandeis University course RSEG176 Cloud Computing.

To Start the Flask server:

1. From a command prompt navigate to the directory that has the file application.py

2. Type the command "python application.py"

3. There will be a response that the server is running on http://127.0.0.1:5000

4. Point a browser to http://127.0.0.1:5000 to run the application

5. Type the command CTRL + c to stop the server


Notes on AWS Elastic Beanstalk Deployment

* Do not need to include .pyc compiled files, just .py

* Make sure requirements.txt file is updated and correct

* Application file must be called application.py

* Deployment upload must be a zip file without the top level directory

* Deployment upload must include .ebextensions directory with a .config file