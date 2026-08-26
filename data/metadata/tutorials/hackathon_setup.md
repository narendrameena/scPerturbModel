
# 🐳 Minimal Jupyter via SSH + Docker (Lambda Cloud)

This guide sets up Jupyter running **inside a Docker container** on your 
Lambda Cloud instance and forwards it securely to your local machine

---

## ✅ 1. SSH into your Lambda instance

On your local machine:

```bash
ssh -i <YOUR_SSH_KEY_PATH> ubuntu@<INSTANCE_IP>
```

Replace `<YOUR_SSH_KEY_PATH>` with your private key path and `<INSTANCE_IP>` with your instance's IP address.

---

## ✅ 2. Start your Docker container with Jupyter

Please replace the image URL with an image of your choosing. The BioNemo image has pytorch and datasets pre-installed 
as well as implementations of several Bio-FMs.
```bash
sudo docker run --gpus all --shm-size=64g -dit \
  --name bionemo \
  -p 8888:8888 \
  -v /home/ubuntu/bionemo_workspace:/workspace \
  nvcr.io/nvidia/clara/bionemo-framework:nightly
```

Enter the container:

```bash
sudo docker exec -it bionemo bash
```

---

## ✅ 3. Install Jupyter and other Python tools inside the container

Inside the container shell, replace the pip and git commands with the packages that you need.

```bash
pip install jupyter anndata
git clone https://huggingface.co/datasets/tahoebio/Tahoe-100M
```

Then launch Jupyter:

```bash
jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --allow-root
```

Note the access token printed in the terminal output (you’ll need it to log in).

---

## ✅ 4. Forward Jupyter port from remote to local

On your **local machine**, open a new terminal and run:
Please note that you can do the same port-forwarding setup to use VScode if you prefer.


```bash
ssh -i <YOUR_SSH_KEY_PATH> -L 8888:127.0.0.1:8888 ubuntu@<INSTANCE_IP>
```

Now visit:

```
http://localhost:8888/?token=<TOKEN>
```

Paste in the token you copied from the container output.

---

## ✅ Optional cleanup

To stop and remove the container:

```bash
sudo docker stop bionemo
sudo docker rm bionemo
```

