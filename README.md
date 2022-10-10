# Redox Web App

[![CI](https://github.com/globus-labs/redox-ai-webapp/actions/workflows/CI.yml/badge.svg)](https://github.com/globus-labs/redox-ai-webapp/actions/workflows/CI.yml)

A web app based on FastAPI that helps us design molecules using AI.

## Installing

Create the environment using Anaconda

```bash
conda env create --file environment.yml
```

## Running

Launch the application using uvicorn:

```commandline
uvicorn redoxweb.app:app --reload
```
