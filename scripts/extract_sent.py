from transformers import AutoTokenizer, AutoModelForSequenceClassification
import os
import pandas as pd
import torch
import preprocessor as p
from tqdm import tqdm

os.chdir("../")
path = os.path.join(os.getcwd(), 'data', 'Bitcoin_tweets.csv')
data = pd.read_csv(path)

data.drop(columns=['user_name', 'user_location', 'user_description', 'user_created',
       'user_followers', 'user_friends', 'user_favourites', 'user_verified',
       'hashtags', 'source', 'is_retweet'], inplace=True)
text = data['text'].fillna("#").tolist()


tokenizer = AutoTokenizer.from_pretrained('nlptown/bert-base-multilingual-uncased-sentiment')

model = AutoModelForSequenceClassification.from_pretrained('nlptown/bert-base-multilingual-uncased-sentiment')

def sentiment_score(review):
    tokens = tokenizer.encode(review, return_tensors='pt')
    result = model(tokens)
    return int(torch.argmax(result.logits))+1

sent=[]
for i in tqdm(text):
    sent.append(sentiment_score(p.clean(i)))

data['sent'] = pd.DataFrame(sent)

data.to_csv('sent.csv')